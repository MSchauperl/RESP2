"""
Unit and regression test for the resp2 package.
"""

# Import package, test suite, and other packages as needed
import resp2
import pytest
import sys

def test_resp2_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "resp2" in sys.modules
