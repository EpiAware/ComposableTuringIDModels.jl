# Issues to file against EpiAware/EpiAwarePackageTools.jl

These are template/scaffold gaps encountered while adopting EpiAwarePackageTools
in EpiAwarePrototype.jl. Each is a candidate GitHub issue. The minimal local
workaround actually applied is noted under each.

---

## 1. `scaffold` writes an MIT `LICENSE`; no way to request a different licence

**Expected:** A package that must ship under Apache-2.0 (e.g. because it
incorporates Apache-2.0 code) can scaffold without having its `LICENSE`
overwritten with MIT, or can tell `scaffold`/`update` which licence to write.

**What happened:** `LICENSE` is a MANAGED template
(`Template("LICENSE", "LICENSE", true, true)` in `src/scaffold.jl`) hardcoded to
the MIT text in `templates/LICENSE`. `scaffold` always (re)writes it, and there
is no `license` keyword in `scaffold_inputs`. Worse, because `LICENSE` is
*managed*, the scheduled `update(pkgdir)` template-sync will silently revert any
package that replaces it with a different licence — re-introducing a licence
incompatibility on every sync.

**Minimal repro:**
```julia
using EpiAwarePackageTools
scaffold("/path/to/pkg")          # writes templates/LICENSE (MIT)
# replace LICENSE with Apache-2.0 by hand
update("/path/to/pkg")            # reverts it back to MIT
```

**Suggested fix:** add a `license` input to `scaffold_inputs` (default `"MIT"`)
selecting among bundled `templates/LICENSE.<spdx>` files, and/or treat `LICENSE`
as package-owned (write-once) rather than managed so a deliberate licence choice
is never reverted by a sync.

**Local workaround applied:** overwrote `LICENSE` with the upstream Apache-2.0
text after scaffolding, added a `NOTICE`, and recorded the override here. The
managed-file revert risk above remains until this is fixed upstream — anyone
running `update()` on this package must re-apply the Apache-2.0 `LICENSE`.
