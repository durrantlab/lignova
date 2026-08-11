# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

# pylint: skip-file

import os
import shutil


def pytest_sessionstart(session):  # pytest_configure(config)
    r"""Called after the Session object has been created and
    before performing collection and entering the run test loop.
    """
    # Creates a tmp directory for writing files.
    test_tmp_path = "./tests/tmp"
    if os.path.exists(test_tmp_path):
        shutil.rmtree(test_tmp_path)
    os.makedirs(test_tmp_path, exist_ok=False)
