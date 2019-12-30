from re import match


def scrub(string):
    """Scrubs user-provided string to prevent database injection.

    Args:
        string: User provided database input.

    Returns:
        The string scrubbed of any trailing spaces, punctuation, or additional text.
    """
    return match(r"([A-Za-z0-9+_-]+)", string.strip()).group(1)
