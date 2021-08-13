class XY:
    def __init__(self, x, y):
        self.x = x
        self.y = y


class XYZ(XY):
    def __init__(self, x, y, z):
        super().__init__(x, y)
        self.z = z


class Obstacle:
    def __init__(self, height=None, array_of_surfaces=None, alpha=None):
        if array_of_surfaces is None:
            array_of_surfaces = []
        self.type = 'screen' if len(array_of_surfaces) == 1 else 'obstacle' \
            if len(array_of_surfaces) == 4 else 'unknown'
        self.height = height
        self.surfaces = array_of_surfaces
        self.alpha = alpha


class Surface:
    def __init__(self, point_a, point_b):
        self.starting_point = point_a
        self.ending_point = point_b
        mid_x = (point_a.x + point_b.x) / 2
        mid_y = (point_a.y + point_b.y) / 2
        self.center = XY(mid_x, mid_y)
        self.a = point_b.y - point_a.y
        self.b = point_a.x - point_b.x
        self.c = (self.a * point_a.x) + (self.b * point_a.y)


def get_distance_between_points(point_a, point_b):
    import math as m

    x2 = (point_b.x - point_a.x) ** 2
    y2 = (point_b.y - point_a.y) ** 2

    if hasattr(point_a, 'z') and hasattr(point_b, 'z'):
        z2 = (point_b.z - point_a.z) ** 2
        return m.sqrt(x2 + y2 + z2)

    return m.sqrt(x2 + y2)


class Source:
    def __init__(self, coords, order, receiver, alphas=None):
        from functools import reduce
        self.coords = coords
        self.order = order
        if alphas is not None:
            alphas_inverse = map(lambda a: 1 - a, alphas)
            self.alpha_factor = reduce(lambda a1, a2: a1 * a2, alphas_inverse)
        else:
            self.alpha_factor = 1
        if alphas is None:
            alphas = []
        self.alphas = alphas
        self.distance_from_receiver = get_distance_between_points(receiver, coords)


def flatten(collection):
    return [item for sublist in collection for item in sublist]


def reflect_point_over_surface(point, surface):
    if surface.a == 0:
        y = surface.c / surface.b
        distance_from_y = point.y - y
        return XYZ(point.x, y - distance_from_y, point.z)
    if surface.b == 0:
        x = surface.c / surface.a
        distance_from_x = point.x - x
        return XYZ(x - distance_from_x, point.y, point.z)
    surface_slope = (- surface.a) / surface.b
    surface_y_intercept = surface.c / surface.b
    perpendicular_slope = -1 / surface_slope
    perpendicular_y_intercept = point.y - (perpendicular_slope * point.x)

    cross_point_x = (perpendicular_y_intercept - surface_y_intercept) / (surface_slope - perpendicular_slope)
    cross_point_y = perpendicular_slope * cross_point_x + perpendicular_y_intercept

    x_diff = point.x - cross_point_x
    y_diff = point.y - cross_point_y

    return XYZ(cross_point_x - x_diff, cross_point_y - y_diff, point.z)


def get_visible_surfaces(source, obstacles):
    import math as m

    def sort_by_first_element(tup):
        return tup[0]

    source_point = XY(source.coords.x, source.coords.y)
    visible_surfaces = []
    for obs in obstacles:
        if obs.type == 'screen':
            visible_surfaces.append((obs.surfaces[0], obs.alpha))
            continue
        distance_from_center = []
        for surf in obs.surfaces:
            distance_from_center.append(
                (m.sqrt((source_point.x - surf.center.x)**2 + (source_point.y - surf.center.y)**2), surf)
            )
        distance_from_center.sort(key=sort_by_first_element)
        visible_surfaces.append((distance_from_center[0][1], obs.alpha))
        visible_surfaces.append((distance_from_center[1][1], obs.alpha))
    return visible_surfaces


def image_source_mtd():
    import math as m
    import matplotlib.pyplot as plt

    receiver = XYZ(2.5, 47.5, 4)

    source = Source(XYZ(42.5, 7.5, 4.5), 0, receiver)
    source_power_lvl = 120
    source_power = 10 ** (source_power_lvl / 10) * (10 ** (-12))

    # max number of reflections
    N = 2
    m_alfa = 0.005066

    obstacles = [
        Obstacle(6, [
            Surface(XY(27.5, 50), XY(42.5, 50)),
            Surface(XY(27.5, 57.5), XY(42.5, 57.5)),
            Surface(XY(27.5, 50), XY(27.5, 57.5)),
            Surface(XY(42.5, 57.5), XY(42.5, 50)),
        ], 0.2),
        Obstacle(8, [
            Surface(XY(5, 22.5), XY(12.5, 15)),
            Surface(XY(12.5, 15), XY(22.5, 25)),
            Surface(XY(22.5, 25), XY(15, 37.5)),
            Surface(XY(15, 37.5), XY(5, 22.5)),
        ], 0.3),
        Obstacle(7, [
            Surface(XY(67.5, 5), XY(75, 45)),
        ], 0.1),
    ]

    # list of lists of sources grouped by order
    # length of sources will be equal to N + 1
    sources = [[source]]

    for n in range(N):
        source_order = n + 1
        image_sources = []
        previous_order_sources = sources[n]
        for s in previous_order_sources:
            visible_surfaces = get_visible_surfaces(s, obstacles)
            for surf, alpha in visible_surfaces:
                image_source_position = reflect_point_over_surface(s.coords, surf)
                image_source = Source(image_source_position, source_order, receiver, s.alphas + [alpha])
                image_sources.append(image_source)
        # extend sources by appending a list of image sources with higher order
        sources = sources + [image_sources]

    flat_sources = flatten(sources)

    intensity_ref = [source.alpha_factor
                     * m.exp(-m_alfa * source.distance_from_receiver)
                     * source_power / (4 * m.pi * source.distance_from_receiver ** 2)
                     for source in flat_sources]
    R = 429
    t = [source.distance_from_receiver / 340.3 for source in flat_sources]
    pressure_ref = [(i * 2 * R) ** (1 / 2) for i in intensity_ref]
    pressure_ref = [pr / (2 * 10 ** (-5)) for pr in pressure_ref]
    spl_ref = [20 * m.log10(pr) for pr in pressure_ref]

    # plot the results
    plt.stem(t, spl_ref)
    plt.xlabel("Czas [s]")
    plt.ylabel("Poziom ciśnienia akustycznego [dB]")
    plt.show()



if __name__ == '__main__':
    image_source_mtd()
