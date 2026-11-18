(function () {
    "use strict";

    // Wide tables scroll inside their own box instead of the page.
    document.querySelectorAll("article table").forEach(function (t) {
        var wrap = document.createElement("div");
        wrap.className = "table-wrap";
        t.parentNode.insertBefore(wrap, t);
        wrap.appendChild(t);
    });

    // Sidebar toggle on narrow screens.
    var toggle = document.querySelector(".menu-toggle");
    if (toggle) {
        toggle.addEventListener("click", function () {
            var open = document.body.classList.toggle("nav-open");
            toggle.setAttribute("aria-expanded", open ? "true" : "false");
        });
    }

    // Bring the current page into view in a long sidebar.
    var current = document.querySelector(".side-nav li.current > a");
    if (current && current.scrollIntoView) {
        current.scrollIntoView({ block: "center" });
    }

    // Client-side search over the prebuilt page index.
    var input = document.getElementById("search");
    var box = document.getElementById("search-results");
    var navRoot = document.querySelector(".side-nav");
    var root = navRoot ? navRoot.getAttribute("data-root") || "" : "";
    var index = window.MV_SEARCH || [];
    if (!input || !box) { return; }

    function snippet(text, pos, len) {
        var start = Math.max(0, pos - 40);
        var end = Math.min(text.length, pos + len + 60);
        return (start > 0 ? "… " : "") + text.slice(start, end) +
            (end < text.length ? " …" : "");
    }

    function run() {
        var q = input.value.trim().toLowerCase();
        box.textContent = "";
        if (q.length < 2) { box.hidden = true; return; }
        var words = q.split(/\s+/);
        var hits = [];
        index.forEach(function (page) {
            var title = page.t.toLowerCase();
            var text = page.x.toLowerCase();
            var score = 0;
            var first = -1;
            for (var i = 0; i < words.length; i++) {
                var inTitle = title.indexOf(words[i]);
                var inText = text.indexOf(words[i]);
                if (inTitle < 0 && inText < 0) { return; }
                if (inTitle >= 0) { score += 10; }
                if (inText >= 0) {
                    score += 1;
                    if (first < 0) { first = inText; }
                }
            }
            hits.push({ page: page, score: score, pos: first });
        });
        hits.sort(function (a, b) { return b.score - a.score; });
        if (!hits.length) {
            var none = document.createElement("p");
            none.textContent = "No pages match.";
            box.appendChild(none);
        }
        hits.slice(0, 8).forEach(function (h) {
            var a = document.createElement("a");
            a.href = root + h.page.u;
            var b = document.createElement("b");
            b.textContent = h.page.t;
            a.appendChild(b);
            if (h.pos >= 0) {
                var s = document.createElement("span");
                s.textContent = snippet(h.page.x, h.pos, words[0].length);
                a.appendChild(s);
            }
            box.appendChild(a);
        });
        box.hidden = false;
    }
    input.addEventListener("input", run);
})();
