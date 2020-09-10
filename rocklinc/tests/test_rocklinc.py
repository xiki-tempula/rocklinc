"""
Unit and regression test for the rocklinc package.
"""

# Import package, test suite, and other packages as needed
import rocklinc
import pytest
import sys

def test_rocklinc_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "rocklinc" in sys.modules
