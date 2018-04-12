import threading
import logging
logging.basicConfig(format='%(asctime)s %(message)s')

from caproto.threading.pyepics_compat import get_pv, caput, caget

thread_class = threading.Thread
pv_form = 'time'


def setup(logger):
    ...
