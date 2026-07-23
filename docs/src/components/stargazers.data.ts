// MANAGED by EpiAwarePackageTools.scaffold — do not edit by hand.
// Build-time loader for the navbar StarUs widget: fetches the repo's live
// GitHub star count. Runs on the docs builder (VitePress data loader), not in
// the browser. A `GITHUB_TOKEN` in the environment raises the API rate limit.
//
// A decorative star count must never fail a docs build. A transient GitHub API
// error (5xx / rate limit / network) is retried, then degraded to an unknown
// count so the page still renders. Only a genuine 404 (misconfigured REPO)
// throws, since re-running will not fix it.
const REPO = "EpiAware/ComposableTuringIDModels.jl";
const MAX_ATTEMPTS = 3;

export default {
  async load() {
    try {
      const {stargazers_count} = await github(`/repos/${REPO}`);
      return stargazers_count;
    } catch (error) {
      if (error.status === 404) throw error;
      console.warn(
        `StarUs: could not fetch the star count (${error.message}); ` +
          "rendering without it."
      );
      return NaN;
    }
  }
};

async function github(
  path,
  {
    authorization = process.env.GITHUB_TOKEN && `token ${process.env.GITHUB_TOKEN}`,
    accept = "application/vnd.github.v3+json"
  } = {}
) {
  const url = new URL(path, "https://api.github.com");
  const headers = {...(authorization && {authorization}), accept};
  for (let attempt = 1; ; attempt++) {
    let response;
    try {
      response = await fetch(url, {headers});
    } catch (error) {
      // Network-level failure: retry, then let the caller degrade.
      if (attempt >= MAX_ATTEMPTS) throw error;
      await sleep(500 * 2 ** (attempt - 1));
      continue;
    }
    if (response.ok) return await response.json();
    // Retry transient server / rate-limit errors; surface the rest.
    const transient = response.status === 429 || response.status >= 500;
    if (!transient || attempt >= MAX_ATTEMPTS) {
      const error = new Error(`fetch error: ${response.status} ${url}`);
      error.status = response.status;
      throw error;
    }
    await sleep(500 * 2 ** (attempt - 1));
  }
}

function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}
