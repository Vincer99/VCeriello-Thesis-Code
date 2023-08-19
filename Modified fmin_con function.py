# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 11:22:28 2023

@author: vince
"""

def fmin_con(objective_function, x0, sigma0,
             g=no_constraints, g2 = no_constraints, g3=no_constraints, h=no_constraints, post_optimization=False,
             archiving=True, **kwargs):

    if 'parallel_objective' in kwargs:
        raise ValueError("`parallel_objective` parameter is not supported by cma.fmin_con")
    if post_optimization and h != no_constraints and (
            not isinstance(post_optimization, float) or post_optimization <= 0):
        raise ValueError("When equality constraints are given, the argument"
                         "``post_optimization`` must be a strictly positive "
                         "float indicating the error on the inequality constraints")
    # prepare callback list
    if callable(kwargs.setdefault('callback', [])):
        kwargs['callback'] = [kwargs['callback']]

    global _al  # for debugging, may be removed at some point
    F = []
    G = []
    _al = AugmentedLagrangian(len(x0))
    _al_set_logging(_al, kwargs)

    # _al.chi_domega = 1.1
    # _al.dgamma = 1.5

    best_feasible_solution = ot.BestSolution2()
    if archiving:
        archives = [
            _constraints_handler.ConstrainedSolutionsArchive(_constraints_handler._g_pos_max),
            _constraints_handler.ConstrainedSolutionsArchive(_constraints_handler._g_pos_sum),
            _constraints_handler.ConstrainedSolutionsArchive(_constraints_handler._g_pos_squared_sum),
        ]
    else:
        archives = []

    def f(x):
        F.append(objective_function(x))
        return F[-1]
    def constraints(x):
        gvals, g2vals, g3vals, hvals = g(x), g2(x), g3(x), h(x)
        # set m and equality attributes of al
        if _al.lam is None:  # TODO: better abide by an "official" interface?
            _al.set_m(len(gvals) + len(g2vals) + len(g3vals) + len(hvals))
            _al._equality = np.asarray(len(gvals) * [False] + len(g2vals) * [False] + len(g3vals) * [False] + len(hvals) * [True],
                                       dtype='bool')
        G.append(list(gvals) + list(g2vals)  + list(g3vals) + list(hvals))
        return G[-1]
    def auglag(x):
        fval, gvals, g2vals, g3vals = f(x), constraints(x), constraints(x), constraints(x)
        alvals = _al(gvals)
        al2vals = _al(g2vals)
        al3vals=_al(g3vals)
        if all([gi <= 0 for gi in gvals]):
            best_feasible_solution.update(fval, x,
                info={'x':x, 'f': fval, 'g':gvals, 'g_al':alvals})
        if all([gi <= 0 for gi in g2vals]):
            best_feasible_solution.update(fval, x,
                info={'x':x, 'f': fval, 'g2':g2vals, 'g2_al':al2vals})
        if all([gi <= 0 for gi in g3vals]):
            best_feasible_solution.update(fval, x,
                info={'x':x, 'f': fval, 'g3':g3vals, 'g3_al':al3vals})
        info = _constraints_handler.constraints_info_dict(
                    _al.count_calls, x, fval, gvals, alvals)
        info2 = _constraints_handler.constraints_info_dict(
                    _al.count_calls, x, fval, g2vals, al2vals)
        info3 = _constraints_handler.constraints_info_dict(
                    _al.count_calls, x, fval, g3vals, al3vals)
        for a in archives:
            a.update(fval, gvals, info)
            a.update(fval, g2vals, info2)
            a.update(fval, g3vals, info3)
        return fval + sum(alvals)+ sum(al2vals)+ sum(al3vals)
    def set_coefficients(es):
        _al.set_coefficients(F, G)
        F[:], G[:] = [], []
    def update(es):
        x = es.ask(1, sigma_fac=0)[0]
        _al.update(f(x), constraints(x))

    kwargs['callback'].extend([set_coefficients, update])
    # The smallest observed f-values may be below the limit value f(x^*_feas)
    # because f-values depend on the adaptive multipliers. Hence we overwrite
    # the default tolstagnation value:
    kwargs.setdefault('options', {}).setdefault('tolstagnation', 0)
    _, es = fmin2(auglag, x0, sigma0, **kwargs)
    es.objective_function_complements = [_al]  # for historical reasons only
    es.augmented_lagrangian = _al
    es.best_feasible = best_feasible_solution

    if post_optimization:
        def f_post(x):
            return sum(gi ** 2 for gi in g(x) if gi > 0) + sum(
                       hi * 2 for hi in h(x) if hi * 2 > post_optimization ** 2)
            
        kwargs_post = kwargs.copy()
        kwargs_post.setdefault('options', {})['ftarget'] = 0

        _, es_post = fmin2(f_post, es.result.xfavorite, es.sigma,
                           **kwargs_post)
        if es_post.best.f == 0:
            f = objective_function(es_post.best.x)
            es.best_feasible.update(f, x=es_post.best.x, info={
                'x': es_post.best.x,
                'f': f,
                'g': None  # it's a feasible solution, so we don't really care
            })
            return es.best_feasible.x, es
        x_post = es_post.result.xfavorite
        g_x_post, h_x_post = g(x_post), h(x_post)
        if all([gi <= 0 for gi in g_x_post]) and \
                all([hi * 2 <= post_optimization * 2 for hi in h_x_post]):
            f_x_post = objective_function(x_post)
            es.best_feasible.update(f_x_post, x=x_post, info={
                'x': x_post,
                'f': f_x_post,
                'g': list(g_x_post) + list(h_x_post)
            })
            return x_post, es
        else:
            utils.print_warning('Post optimization was unsuccessful',
                                verbose=es.opts['verbose'])

    es.con_archives = archives
    return es.result.xfavorite, es  # do not return es.best_feasible.x because it could be quite bad
