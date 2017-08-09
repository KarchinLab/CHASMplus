# logging imports
import warnings
import logging
import datetime
import os
import sys

# config handling
import yaml
data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data'))
config_path = os.path.join(data_dir, "config.yaml")
with open(config_path) as handle:
    config_dict = yaml.load(handle)


def start_logging(log_file='', log_level='INFO', verbose=False):
    """Start logging information into the log directory.
    If os.devnull is specified as the log_file then the log file will
    not actually be written to a file.
    """
    if not log_file:
        # create log directory if it doesn't exist
        log_dir = os.path.abspath('log') + '/'
        if not os.path.isdir(log_dir):
            os.mkdir(log_dir)

        # path to new log file
        log_file = log_dir + 'log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'

    # logger options
    lvl = logging.DEBUG if log_level.upper() == 'DEBUG' else logging.INFO

    # ignore warnings if not in debug
    if log_level.upper() != 'DEBUG':
        warnings.filterwarnings('ignore')

    # define logging format
    if verbose:
        myformat = '%(asctime)s - %(name)s - %(levelname)s \n>>>  %(message)s'
    else:
        myformat = '%(message)s'

    # create logger
    if not log_file == 'stdout':
        # normal logging to a regular file
        logging.basicConfig(level=lvl,
                            format=myformat,
                            filename=log_file,
                            filemode='w')
    else:
        # logging to stdout
        root = logging.getLogger()
        root.setLevel(lvl)
        stdout_stream = logging.StreamHandler(sys.stdout)
        stdout_stream.setLevel(lvl)
        formatter = logging.Formatter(myformat)
        stdout_stream.setFormatter(formatter)
        root.addHandler(stdout_stream)
        root.propagate = True
