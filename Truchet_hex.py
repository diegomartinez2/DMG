import numpy as np
import random

from truchet import Truchet

class TruchetHexes(Truchet):
    """A class for creating a Truchet tiling of hexagons."""

    def __init__(self, width, height, s, colour):
        super(TruchetHexes, self).__init__(width, height, s)
        self.colour = colour

    @Truchet.defs_decorator
    def svg_styles(self):
        print('.arc {{ stroke: {}; stroke-width: 2px; fill: none; }}'
              .format(self.colour), file=self.fo)
        print('.hex {{ stroke: {}; stroke-width: 2px; fill: none; }}'
              .format('#eeeeee'), file=self.fo)


    def svg_shape(self, r=None, rule=None):
        """A Truchet figure based on interlinking circular arcs."""

        def arc_path(A, B, r):
            """Semicircular arc path from A=(x0,y0) to B=(x1,y1), radius r."""

            print('<path d="M{},{} A{},{} 0 0 1 {} {}" class="arc"/>'.format(
                  *A, r, r, *B), file=self.fo)

        def line_path(A, B):
            """The straight line "under" the curved arcs.

            NB for now, the "rule" is ignored: only a random weave is produced.

            """

            # Unit vector across the hexagon.
            V = B - A
            V /= np.hypot(*V)
            # Adjust the padding according to the line width to leave a gap
            # either side of the arc weaving "above" it.
            pad = 0.08
            g1 = self.s * ((3 - np.sqrt(3))/2 - pad)
            g2 = self.s * ((3 - np.sqrt(3))/2 + pad)
            Q = np.array([A, A + g1*V, A + V * g2, B - V * g2, B - V * g1, B])
            print('<path d="M{} {} L{} {} M{} {} L{} {} M{} {} L{} {}"'
                  ' class="arc"/>'.format(*Q.ravel()), file=self.fo)


        f1, f2 = 3/2, np.sqrt(3)/2
        if not r:
            r = self.s * f1
        for ix in range(self.nx):
            for iy in range(self.ny):
                # The centre of this hexagon.
                x0, y0 = (ix * self.s * f1,
                          iy * self.s * (f2*2) + self.s*f2*(ix % 2) )
                # The mid-points of each side of the hexagon.
                P = np.empty((6,2))
                P[0] = 0, -self.s * f2
                R = np.array(((0.5, -f2), (f2, 0.5)))
                P[0] = (0, - self.s * f2)
                for i in range(1,6):
                    P[i] = R @ P[i-1]
                P += (x0, y0)

                # If we're drawing the hexagons themselves, these are their
                # vertices.
                Q = np.empty((6,2))
                Q[0] = self.s/2, self.s * f2
                for i in range(1,6):
                    Q[i] = R @ Q[i-1]
                Q += (x0, y0)

                # If drawing the hexagons, uncomment these lines.
                #print('<path d="M{} {} L{} {} L{} {} L{} {} L{} {}'
                #       ' L{} {}z" class="hex"/>'.format(*Q.ravel()),
                #       file=self.fo)

                # Randomly orient the hexagon by cyclicly shifting the
                # coordinate rows 0, 1 or 2 times.
                P = np.roll(P, random.randint(0,2), axis=0)
                # Draw the arcs and line.
                arc_path(P[0], P[4], r)
                arc_path(P[3], P[1], r)
                line_path(P[2], P[5])

if __name__ == '__main__':
    truchet = TruchetHexes(800, 800, 25, colour='#4f3e90')
    truchet.make_svg('hexes.svg')
