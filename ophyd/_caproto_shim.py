import threading
import logging
logging.basicConfig(format='%(asctime)s %(message)s')

from caproto.threading.client import PVContext

ctx = PVContext(log_level='DEBUG')

def get_pv(*args, **kwargs):
    return ctx.get_pv(*args, **kwargs)


def caget(pvname, *args, **kwargs):
    return get_pv(pvname).get(*args, **kwargs)


def caput(pvname, *args, **kwargs):
    return get_pv(pvname).put(*args, **kwargs)

thread_class = threading.Thread
pv_form = 'time'


def setup(logger):
    ...
