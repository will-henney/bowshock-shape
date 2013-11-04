import datetime
import getpass
import sys
import os

def run_info():
    """
    Return info about this run:
    1. Name of program
    2. User who is running it
    3. Date and time it was run
    """
    return "{}: <{}> {}".format(
        os.path.basename(sys.argv[0]), 
        getpass.getuser(), 
        datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
    )
