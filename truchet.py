class Truchet:
    """Base class for a Truchet tiling."""

    def __init__(self, width, height, s):
        """Initialize the class with image size and tile size, s."""

        self.width, self.height = width, height
        self.s = s
        self.nx, self.ny = int(width // s), int(height // s)
        self.fo = None

    def preamble(self):
        """The usual SVG preamble, including the image size."""

        print('<?xml version="1.0" encoding="utf-8"?>\n'

        '<svg xmlns="http://www.w3.org/2000/svg"\n' + ' '*5 +
           'xmlns:xlink="http://www.w3.org/1999/xlink" width="{}" height="{}" >'
                .format(self.width, self.height), file=self.fo)

    def defs_decorator(func):
        """For convenience, wrap the CSS styles with the needed SVG tags."""

        def wrapper(self):
            print("""
            <defs>
            <style type="text/css"><![CDATA[""", file=self.fo)

            func(self)

            print("""]]></style>
            </defs>""", file=self.fo)
        return wrapper

    def svg_shape(self, *args, **kwargs):
        """Override this function in the derived class."""

    def make_svg(self, filename, *args, **kwargs):
        """Create the tiling image as an SVG file with name filename.

        Custom arguments are passed to the derived class's svg_shape method.

        """

        self.fo = open(filename, 'w')
        self.preamble()
        self.svg_styles()
        self.svg_shape(*args, **kwargs)
        print('</svg>', file=self.fo)
