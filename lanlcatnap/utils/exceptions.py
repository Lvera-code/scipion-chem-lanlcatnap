"""Exception hierarchy for this plugin: never let a raw
FileNotFoundError/ValueError escape to the Scipion GUI without an
actionable message.
"""


class LANLCATNAPExecutionError(Exception):
    """Failed to cross-reference against the LANL/CATNAP databases: missing
    installation, or the reference CSV/TSV does not match the expected
    format."""
