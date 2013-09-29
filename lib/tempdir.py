import os
import tempfile
import shutil
import contextlib

@contextlib.contextmanager
def temporary_directory():
    """A context manager which changes the working directory to a
    temporary directory, then cleans it up and cds back
    """
    prev_cwd = os.getcwd()
    tempdir = tempfile.mkdtemp()
    os.chdir(tempdir)
    yield
    os.chdir(prev_cwd)
    shutil.rmtree(tempdir)
