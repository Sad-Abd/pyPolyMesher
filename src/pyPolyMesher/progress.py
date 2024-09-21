"""
Progress Bar Module

This module provides a customizable progress bar for command-line interfaces.
It displays completion percentage, elapsed time, and estimated time to completion.

Adapted from https://github.com/fogleman
"""

import sys
import time


def pretty_time(seconds):
    """
    Convert seconds to a formatted string of hours, minutes, and seconds.

    Args:
        seconds (float): The number of seconds to convert.

    Returns:
        str: A string in the format "H:MM:SS".
    """
    seconds = int(round(seconds))
    s = seconds % 60
    m = (seconds // 60) % 60
    h = seconds // 3600
    return "%d:%02d:%02d" % (h, m, s)


class Bar(object):
    """
    A customizable progress bar for command-line interfaces.

    This class creates and updates a progress bar, showing completion percentage,
    current value, elapsed time, and estimated time to completion.
    """

    def __init__(self, max_value=100, min_value=0, enabled=True):
        """
        Initialize the progress bar.

        Args:
            max_value (int): The maximum value of the progress bar. Defaults to 100.
            min_value (int): The minimum value of the progress bar. Defaults to 0.
            enabled (bool): Whether the progress bar should be displayed. Defaults to True.
        """
        self.min_value = min_value
        self.max_value = max_value
        self.value = min_value
        self.start_time = time.time()
        self.enabled = enabled

        self.error = 0
        self.min_error = 100

    @property
    def percent_complete(self):
        """Calculate and return the percentage of completion."""
        t = (self.value - self.min_value) / (self.max_value - self.min_value)
        return t * 100

    @property
    def elapsed_time(self):
        """Calculate and return the elapsed time since the start."""
        return time.time() - self.start_time

    @property
    def eta(self):
        """Estimate and return the time remaining to completion."""
        t = self.percent_complete / 100
        if t == 0:
            return 0
        return (1 - t) * self.elapsed_time / t

    def increment(self, delta, error):
        """
        Increment the current value and update the error.

        Args:
            delta (float): The amount to increment the current value by.
            error (float): The current error value.
        """
        self.update(self.value + delta)
        self.error = error
        if error < self.min_error:
            self.min_error = error

    def update(self, value):
        """
        Update the current value and redraw the progress bar.

        Args:
            value (float): The new current value.
        """
        self.value = value
        if self.enabled:
            sys.stdout.write("  %s    \r" % self.render())
            sys.stdout.flush()

    def done(self):
        """Mark the task as complete and stop the progress bar."""
        self.update(self.max_value)
        self.stop()

    def stop(self):
        """Stop the progress bar and move to a new line."""
        if self.enabled:
            sys.stdout.write("\n")
            sys.stdout.flush()

    def render(self):
        """Render the complete progress bar as a string."""
        items = [
            self.render_percent_complete(),
            self.render_value(),
            self.render_bar(),
            self.render_elapsed_time(),
            self.render_eta(),
            self.render_error(),
        ]
        return " ".join(items)

    def render_percent_complete(self):
        """Render the percentage complete as a string."""
        return "%3.0f%%" % self.percent_complete

    def render_value(self):
        """Render the current value as a string."""
        if self.min_value == 0:
            return "(%g of %g)" % (self.value, self.max_value)
        else:
            return "(%g)" % (self.value)

    def render_bar(self, size=30):
        """
        Render the visual progress bar.

        Args:
            size (int): The total width of the progress bar. Defaults to 30.

        Returns:
            str: A string representation of the progress bar.
        """
        a = int(round(self.percent_complete / 100.0 * size))
        b = size - a
        return "[" + "#" * a + "-" * b + "]"

    def render_elapsed_time(self):
        """Render the elapsed time as a formatted string."""
        return pretty_time(self.elapsed_time)

    def render_eta(self):
        """Render the estimated time to completion as a formatted string."""
        return pretty_time(self.eta)

    def render_error(self):
        """Render the current and minimum error as a string."""
        return "[Error: %f Min Error: %f]" % (self.error, self.min_error)
