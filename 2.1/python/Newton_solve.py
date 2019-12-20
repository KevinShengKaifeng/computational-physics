def D(x, func, delta=1e-3):
    return (func(x+delta/2)-func(x-delta/2))/delta


def newton_solve_fixed_x0(func, x0):
    max_step = 200
    i = 0
    deltamin = 1e-8
    while i < max_step:
        d = D(x0, func)
        if abs(d) < 1e-3:
            d = 5e-2
        deltax = -func(x0) / d
        x0 += deltax
        i += 1
        if abs(deltax) < deltamin:
            return x0
    return


def newton_solve(func, x0=0, root_range=()):
    solve = newton_solve_fixed_x0(func, x0)
    if solve is not None:
        return solve
    elif not root_range:
        print("root not found!")
        return
    find_time = 10
    for i in range(find_time):
        solve = newton_solve_fixed_x0(func, root_range[0]+(root_range[1]-root_range[0])/find_time*i)
        if solve is not None:
            return solve
    print("root not found!")
    return