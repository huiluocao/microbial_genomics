#!/usr/bin/env python3

import os
import sys
import base64
import argparse
import logging
import urllib.request
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode


SERVER_ADDRESS = "https://enterobase.warwick.ac.uk"
DATABASE = "clostridium"


def create_request(url, api_token):
    auth_string = f"{api_token}: "
    base64string = base64.b64encode(auth_string.encode("utf-8")).decode("utf-8")

    headers = {
        "Authorization": f"Basic {base64string}"
    }

    return urllib.request.Request(url, None, headers)


def read_barcodes_from_file(barcode_file):
    barcodes = []

    with open(barcode_file, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            barcodes.append(line)

    return barcodes


def download_assembly(assembly_barcode, api_token, outdir):
    params = {
        "assembly_barcode": assembly_barcode,
        "database": DATABASE
    }

    url = f"{SERVER_ADDRESS}/upload/download?{urlencode(params)}"

    outfile = os.path.join(outdir, f"{assembly_barcode}.fasta")

    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        logging.info(f"Already exists, skipping: {outfile}")
        return True

    logging.info(f"Downloading {assembly_barcode}")
    logging.info(f"URL: {url}")

    try:
        response = urlopen(create_request(url, api_token))

        content = response.read()

        if not content:
            logging.error(f"No content downloaded for {assembly_barcode}")
            return False

        with open(outfile, "wb") as f:
            f.write(content)

        logging.info(f"Saved: {outfile}")
        return True

    except HTTPError as e:
        error_text = e.read().decode(errors="ignore")
        logging.error(
            f"HTTP error while downloading {assembly_barcode}: "
            f"{e.code} {e.reason} <{e.geturl()}> {error_text}"
        )
        return False

    except URLError as e:
        logging.error(f"URL error while downloading {assembly_barcode}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Download assemblies from EnteroBase clostridium database using assembly barcodes."
    )

    parser.add_argument(
        "-b",
        "--barcode",
        nargs="*",
        help="One or more EnteroBase assembly barcodes, e.g. CLO_AA0007AA_AS"
    )

    parser.add_argument(
        "-f",
        "--barcode-file",
        help="Text file with one assembly barcode per line."
    )

    parser.add_argument(
        "-o",
        "--outdir",
        default="assemblies",
        help="Output directory. Default: assemblies"
    )

    parser.add_argument(
        "--token",
        default=os.getenv("ENTEROBASE_API_TOKEN"),
        help="EnteroBase API token. Prefer using ENTEROBASE_API_TOKEN environment variable."
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s"
    )

    if not args.token:
        logging.error("No API token provided. Set ENTEROBASE_API_TOKEN or use --token.")
        sys.exit(1)

    barcodes = []

    if args.barcode:
        barcodes.extend(args.barcode)

    if args.barcode_file:
        barcodes.extend(read_barcodes_from_file(args.barcode_file))

    barcodes = list(dict.fromkeys(barcodes))

    if not barcodes:
        logging.error("No assembly barcodes provided.")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    success = 0
    failed = 0

    for barcode in barcodes:
        ok = download_assembly(
            assembly_barcode=barcode,
            api_token=args.token,
            outdir=args.outdir
        )

        if ok:
            success += 1
        else:
            failed += 1

    logging.info("Download complete.")
    logging.info(f"Successful downloads: {success}")
    logging.info(f"Failed downloads: {failed}")


if __name__ == "__main__":
    main()
