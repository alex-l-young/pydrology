# =======================================================
# NEXRAD Level-II Scan Requests
# =======================================================

# Library imports.
from pathlib import Path
import nexradaws
import six
import threading
from queue import Queue
import logging
from time import time
import argparse

# Logging.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

logger = logging.getLogger(__name__)

# =======================================================
# NEXRAD Doppler Data.
# =======================================================
def download_nexrad_scan(scan, output_file):
    conn = nexradaws.NexradAwsInterface()
    localfiles = conn.download(scan, output_file)
    six.print_(localfiles.success)
    six.print_(localfiles.success[0].filepath)

    return localfiles


def setup_download_dir(download_dir):
    if not download_dir.exists():
        download_dir.mkdir()
    return download_dir


def get_scans(date, site):
    # Split the date into year, month, and day.
    date_split = date.split('-')
    year, month, day = (date_split[0], date_split[1], date_split[2])

    conn = nexradaws.NexradAwsInterface()
    scans = conn.get_avail_scans(year, month, day, site)

    scans_to_download = []
    for scan in scans:
        scan_str = str(scan)
        scan_split = scan_str.split('_')

        # Don't process MDM files.
        if 'MDM' not in scan_split[-1]:
            scans_to_download.append(scan)
        else:
            continue

    print(f'Found {len(scans_to_download)} scans for this date...', '\n')

    return scans_to_download


class DownloadWorker(threading.Thread):

    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.queue = queue

    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            scan, output_file = self.queue.get()
            try:
                download_nexrad_scan(scan, output_file)
            finally:
                self.queue.task_done()


def main(output_dir, date, site, threads):
    ts = time()
    download_dir = output_dir / date
    download_dir = setup_download_dir(download_dir)
    scans = get_scans(date, site)
    # Create a queue to communicate with the worker threads
    queue = Queue()
    # Create 8 worker threads
    for x in range(threads):
        worker = DownloadWorker(queue)
        # Setting daemon to True will let the main thread exit even though the workers are blocking
        worker.daemon = True
        worker.start()
    # Put the tasks into the queue as a tuple
    for scan in scans:
        logger.info('Queueing {}'.format(scan))
        queue.put((scan, download_dir))
    # Causes the main thread to wait for the queue to finish processing all the tasks
    queue.join()
    logging.info('Took %s', time() - ts)


def parse_arguments():
    # Parse the arguments.
    parser = argparse.ArgumentParser(description='Download NEXRAD Scans from a single day.')
    parser.add_argument('--dir', help='Output directory. Files will be saved in sub-directory using the date.')
    parser.add_argument('--threads', help='Number of download threads to use at once.')
    parser.add_argument('--date', help='Date of the scans yyyy-mm-dd')
    parser.add_argument('--site', help='NEXRAD site code.')
    args = parser.parse_args()

    dir = Path(args.dir)
    site = args.site
    date = args.date
    threads = int(args.threads)

    parsed_args = [dir, date, site, threads]

    return parsed_args


if __name__ == '__main__':
    parsed_args = parse_arguments()
    main(*parsed_args)

