window.MathJax = {
  "tex": {
    "packages": [
      "base",
      "ams",
      "autoload",
      "configmacros",
      "physics"
    ],
    "macros": {
      "x": "\\bm{\\text{x}}",
      "jump": [
        "\\left[\\!\\left[ #1 \\right]\\!\\right]",
        1
      ],
      "R": "\\mathbb{R}",
      "sign": "\\operatorname{sign}",
      "s": "\\sigma",
      "phi": "\\varphi",
      "div": [
        "\\operatorname{div}\\left(#1\\right)",
        1
      ],
      "eps": "\\varepsilon",
      "bm": [
        "{\\boldsymbol{#1}}",
        1
      ],
      "esssup": "\\operatorname*{ess\\,sup}",
      "essinf": "\\operatorname*{ess\\,inf}",
      "Cell": "\\mathcal{C}",
    },
    "tags": "ams"
  },
  "loader": {
    "load": [
      "[tex]/physics"
    ]
  },
  startup: {
    ready: () => {
      MathJax.startup.defaultReady();
      MathJax.startup.promise.then(() => {
        subscribe();
      });
    },
  },
};
  

function subscribe() { document$.subscribe(() => { 
  MathJax.startup.output.clearCache()
  MathJax.typesetClear()
  MathJax.texReset()
  MathJax.typesetPromise().then(() => {
    document.dispatchEvent(new Event('MathJaxFinished'));
  });
})}