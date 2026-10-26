using TestItemRunner
pats = ARGS
@run_package_tests filter = ti -> !(:ad in ti.tags) &&
    (isempty(pats) || any(p -> occursin(p, ti.filename), pats))
