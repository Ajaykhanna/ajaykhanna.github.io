// Unified scripts (projects + enhancements). Load-guarded to avoid duplicate execution.
if (window.__unifiedScriptsLoaded) {
  console.warn("enhanced-scripts.js: already loaded — skipping.");
} else {
  window.__unifiedScriptsLoaded = true;

  (function () {
    "use strict";

    /* ---------- Initializers ---------- */
    function initAll() {
      initParticles();
      initTypedText();
      initAOS();
      initCounters();
      initProgressBar();
      initSmoothScroll();
      initLazyLoading();
      initNavToggle();
      initFAB();
      // slider init handled inside their scoped modules
    }

    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", initAll);
    } else {
      initAll();
    }

    /* ---------- Particles ---------- */
    function initParticles() {
      if (typeof particlesJS === "undefined") return;
      try {
        particlesJS("particles-js", {
          particles: {
            number: { value: 80, density: { enable: true, value_area: 800 } },
            color: { value: ["#3b82f6", "#10b981", "#14b8a6", "#8b5cf6"] },
            shape: { type: "circle" },
            opacity: { value: 0.5, random: true },
            size: { value: 3, random: true },
            line_linked: {
              enable: true,
              distance: 150,
              color: "#ffffff",
              opacity: 0.2,
              width: 1,
            },
            move: { enable: true, speed: 2, out_mode: "out" },
          },
          interactivity: {
            detect_on: "canvas",
            events: {
              onhover: { enable: true, mode: "grab" },
              onclick: { enable: true, mode: "push" },
              resize: true,
            },
            modes: { grab: { distance: 140 } },
          },
          retina_detect: true,
        });
      } catch (e) {
        /* ignore */
      }
    }

    /* ---------- Typed.js ---------- */
    function initTypedText() {
      if (typeof Typed === "undefined") return;
      const el = document.getElementById("typed-output");
      if (!el) return;
      new Typed("#typed-output", {
        strings: [
          "Drug Discovery Expert",
          "Energy Transfer Researcher",
          "Battery Materials Scientist",
          "Green Energy Innovator",
          "ML for Chemistry Enthusiast",
          "QM/MM Specialist",
        ],
        typeSpeed: 50,
        backSpeed: 30,
        backDelay: 2000,
        loop: true,
        showCursor: true,
        cursorChar: "|",
      });
    }

    /* ---------- AOS ---------- */
    function initAOS() {
      if (typeof AOS === "undefined") return;
      AOS.init({
        duration: 1000,
        once: false,
        mirror: true,
        offset: 100,
        easing: "ease-in-out",
        anchorPlacement: "top-bottom",
      });
    }

    /* ---------- Counters ---------- */
    function initCounters() {
      const counters = document.querySelectorAll(".counter");
      if (!counters.length) return;
      function animateCounter(counter) {
        const target = parseInt(counter.dataset.target || "0", 10);
        const duration = 2000,
          steps = 60,
          increment = target / steps;
        let current = 0;
        const t = setInterval(() => {
          current += increment;
          if (current >= target) {
            counter.textContent = target + "+";
            clearInterval(t);
          } else counter.textContent = Math.floor(current);
        }, duration / steps);
      }
      const obs = new IntersectionObserver(
        (entries) => {
          entries.forEach((e) => {
            if (e.isIntersecting && !e.target.classList.contains("counted")) {
              e.target.classList.add("counted");
              animateCounter(e.target);
            }
          });
        },
        { threshold: 0.5 }
      );
      counters.forEach((c) => obs.observe(c));
    }

    /* ---------- Progress bar ---------- */
    function initProgressBar() {
      const progressBar = document.getElementById("progress-bar");
      if (!progressBar) return;
      window.addEventListener(
        "scroll",
        () => {
          const h = document.documentElement;
          const percent =
            ((window.pageYOffset || h.scrollTop) /
              (h.scrollHeight - window.innerHeight)) *
            100;
          progressBar.style.width = Math.min(Math.max(percent, 0), 100) + "%";
        },
        { passive: true }
      );
    }

    /* ---------- Smooth scroll ---------- */
    function initSmoothScroll() {
      document.querySelectorAll('a[href^="#"]').forEach((a) => {
        a.addEventListener("click", function (e) {
          const href = this.getAttribute("href");
          if (!href || href === "#") return;
          const target = document.querySelector(href);
          if (!target) return;
          e.preventDefault();
          const navH = document.querySelector("nav")?.offsetHeight || 0;
          const pos =
            target.getBoundingClientRect().top + window.pageYOffset - navH;
          window.scrollTo({ top: pos, behavior: "smooth" });
        });
      });
    }

    /* ---------- Lazy loading ---------- */
    function initLazyLoading() {
      if (!("IntersectionObserver" in window)) return;
      const obs = new IntersectionObserver((entries, o) => {
        entries.forEach((en) => {
          if (!en.isIntersecting) return;
          const img = en.target;
          img.src = img.dataset.src || img.src;
          img.classList.remove("lazy");
          o.unobserve(img);
        });
      });
      document.querySelectorAll("img.lazy").forEach((i) => obs.observe(i));
    }

    /* ---------- Nav toggle (safe) ---------- */
    function initNavToggle() {
      const btn = document.getElementById("nav-toggle");
      const content = document.getElementById("nav-content");
      if (!btn || !content) return;
      btn.addEventListener("click", () => content.classList.toggle("hidden"));
    }

    /* ---------- Sliders (scoped, expose functions) ---------- */
    (function projectSliderModule() {
      if (typeof window.projectSliderPos === "undefined")
        window.projectSliderPos = 0;
      const el = document.getElementById("projectSlider");
      if (!el) return;
      function slideProject(direction) {
        const w = el.children[0]?.offsetWidth || 0;
        const total = el.children.length || 1;
        window.projectSliderPos = Math.max(
          Math.min(window.projectSliderPos + direction, 0),
          -total + 1
        );
        el.style.transform = `translateX(${window.projectSliderPos * w}px)`;
      }
      window.slideProject = slideProject;
      window.addEventListener("resize", () => {
        const w = el.children[0]?.offsetWidth || 0;
        el.style.transform = `translateX(${window.projectSliderPos * w}px)`;
      });
    })();

    (function tutorialSliderModule() {
      if (typeof window.tutorialSliderPos === "undefined")
        window.tutorialSliderPos = 0;
      const el = document.getElementById("tutorialSlider");
      if (!el) return;
      function slideTutorials(direction) {
        const w = el.children[0]?.offsetWidth || 0;
        const total = el.children.length || 1;
        window.tutorialSliderPos = Math.max(
          Math.min(window.tutorialSliderPos + direction, 0),
          -total + 1
        );
        el.style.transform = `translateX(${window.tutorialSliderPos * w}px)`;
      }
      window.slideTutorials = slideTutorials;
      window.addEventListener("resize", () => {
        const w = el.children[0]?.offsetWidth || 0;
        el.style.transform = `translateX(${window.tutorialSliderPos * w}px)`;
      });
    })();

    /* ---------- FAB (scoped) ---------- */
    function initFAB() {
      const fabMain = document.getElementById("fabMain");
      const fabActions = document.getElementById("fabActions");
      const fabSocialMenu = document.getElementById("fabSocialMenu");
      const fabActionButtons = document.querySelectorAll(".fab-action");
      if (!fabMain || !fabActions) return;
      let isMenuOpen = false,
        isSocialOpen = false;
      function closeSocial() {
        isSocialOpen = false;
        fabSocialMenu?.classList.remove("active");
      }
      function closeAll() {
        isMenuOpen = false;
        isSocialOpen = false;
        fabMain.classList.remove("active");
        fabActions.classList.remove("active");
        fabSocialMenu?.classList.remove("active");
      }
      function toggle() {
        isMenuOpen = !isMenuOpen;
        fabMain.classList.toggle("active");
        fabActions.classList.toggle("active");
        if (!isMenuOpen) closeSocial();
      }
      fabMain.addEventListener("click", (e) => {
        e.stopPropagation();
        toggle();
      });
      fabActionButtons.forEach((btn) =>
        btn.addEventListener("click", function (e) {
          e.stopPropagation();
          const action = this.dataset.action;
          if (action === "scroll-top") {
            window.scrollTo({ top: 0, behavior: "smooth" });
            closeAll();
            return;
          }
          if (action === "contact") {
            const sec = document.getElementById("contact");
            if (sec) {
              const navH = document.querySelector("nav")?.offsetHeight || 0;
              const pos =
                sec.getBoundingClientRect().top + window.pageYOffset - navH;
              window.scrollTo({ top: pos, behavior: "smooth" });
            }
            closeAll();
            return;
          }
          if (action === "resume") {
            window.open("./resume.html", "_blank");
            closeAll();
            return;
          }
          if (action === "social") {
            isSocialOpen = !isSocialOpen;
            fabSocialMenu?.classList.toggle("active");
            fabActions.classList.remove("active");
          }
        })
      );
      document.addEventListener("click", (e) => {
        if (!e.target.closest(".fab-container")) closeAll();
      });
      document.addEventListener("keydown", (e) => {
        if (e.key === "Escape") closeAll();
      });
      let lastScroll = 0;
      window.addEventListener(
        "scroll",
        () => {
          const s = window.pageYOffset || document.documentElement.scrollTop;
          if (s > lastScroll && s > 300) {
            fabMain.style.transform = "translateY(100px)";
            closeAll();
          } else fabMain.style.transform = "translateY(0)";
          lastScroll = s;
        },
        { passive: true }
      );
    }

    /* ---------- Update active nav ---------- */
    function updateActiveNav() {
      const sections = document.querySelectorAll("section[id]");
      const navLinks = document.querySelectorAll('nav a[href^="#"]');
      let current = "";
      sections.forEach((s) => {
        if (window.pageYOffset >= s.offsetTop - 200) current = s.id;
      });
      navLinks.forEach((l) =>
        l.classList.toggle("active", l.getAttribute("href") === "#" + current)
      );
    }
    window.addEventListener("scroll", updateActiveNav, { passive: true });

    /* ---------- Publication Timeline Filtering ---------- */
    window.filterPublicationsByYear = function (year) {
      const publications = document.querySelectorAll('.publication-item');
      const buttons = document.querySelectorAll('.timeline-year');

      // Update active button
      buttons.forEach(btn => {
        btn.classList.toggle('active', btn.dataset.year === year);
      });

      // Filter publications
      publications.forEach(pub => {
        if (year === 'all' || pub.dataset.year === year) {
          pub.classList.remove('hidden');
        } else {
          pub.classList.add('hidden');
        }
      });
    };

    /* ---------- Collaboration Map Initialization ---------- */
    function initCollaborationMap() {
      const mapElement = document.getElementById('collaboration-map');
      if (!mapElement || typeof L === 'undefined') return;

      try {
        // Initialize map centered on USA
        const map = L.map('collaboration-map').setView([39.8283, -98.5795], 4);

        // Add tile layer
        L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
          attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
          maxZoom: 18,
        }).addTo(map);

        // Collaboration locations
        const locations = [
          {
            name: 'Los Alamos National Laboratory (LANL)',
            coords: [35.8800, -106.3031],
            city: 'Los Alamos, NM, USA',
            publications: ['Chiroptical Response (2025)', 'Perylene Diimide Trimers (2025)'],
            collaborator: 'Sergei Tretiak',
            isHome: true
          },
          {
            name: 'UC Merced (UCM)',
            coords: [37.3661, -120.4250],
            city: 'Merced, CA, USA',
            publications: ['Optical Spectroscopy (2021)', 'Franck-Condon Methods (2024)'],
            collaborator: 'Christine M. Isborn, Sapana V. Shedge'
          },
          {
            name: 'UC Irvine (UCI)',
            coords: [33.6405, -117.8443],
            city: 'Irvine, CA, USA',
            publications: ['Chiroptical Response (2025)'],
            collaborator: 'Shaul Mukamel, Victor M. Freixas'
          },
          {
            name: 'UC San Diego (UCSD)',
            coords: [32.8801, -117.2340],
            city: 'San Diego, CA, USA',
            publications: ['Molecular Polariton (2022)'],
            collaborator: 'Joel Yuen-Zhou'
          },
          {
            name: 'Oregon State University (OSU)',
            coords: [44.5646, -123.2620],
            city: 'Corvallis, OR, USA',
            publications: ['Optical Spectroscopy (2021)', 'Franck-Condon Methods (2024)'],
            collaborator: 'Tim J. Zuehlsdorff'
          },
          {
            name: 'Penn State University',
            coords: [40.7982, -77.8599],
            city: 'State College, PA, USA',
            publications: ['Molecular Polariton (2022)'],
            collaborator: 'Noel C. Giebink, Chiao-Yu Cheng, Nina Krainova, Alyssa Brigeman'
          },
          {
            name: 'University of Bologna',
            coords: [44.4949, 11.3426],
            city: 'Bologna, Italy',
            publications: ['Chiroptical Response (2025)'],
            collaborator: 'Marco Garavelli'
          },
          {
            name: 'Pacific Northwest National Laboratory (PNNL)',
            coords: [46.3458, -119.2781],
            city: 'Richland, WA, USA',
            publications: ['Chiroptical Response (2025)'],
            collaborator: 'Niranjan Govind, Lei Xu'
          },
          {
            name: 'Argonne National Laboratory (ANL)',
            coords: [41.7123, -87.9850],
            city: 'Lemont, IL, USA',
            publications: ['Chiroptical Response (2025)'],
            collaborator: 'Jérémy R. Rouxel'
          },
          {
            name: 'University of New Mexico (UNM)',
            coords: [35.0844, -106.6504],
            city: 'Albuquerque, NM, USA',
            publications: ['Perylene Diimide Trimers (2025)'],
            collaborator: 'Jean-Hubert Olivier'
          },
          {
            name: 'Universidad Nacional de Quilmes/CONICET',
            coords: [-34.7153, -58.2527],
            city: 'Bernal, Argentina',
            publications: ['Perylene Diimide Trimers (2025)'],
            collaborator: 'Sebastian Fernandez-Alberti'
          },
          {
            name: 'National Institute of Technology Rourkela',
            coords: [22.2500, 84.9000],
            city: 'Odisha, India',
            publications: ['Ligand Driven Electron Counting Rule (2018)'],
            collaborator: 'Rakesh Parida, G. Naaresh Reddy, Gourisankar Roymahapatra, Santanab Giri'
          }
        ];

        // Add markers
        locations.forEach(loc => {
          const iconColor = loc.isHome ? 'red' : 'blue';
          const iconSize = loc.isHome ? [12, 12] : [8, 8];

          const marker = L.circleMarker(loc.coords, {
            radius: loc.isHome ? 10 : 6,
            fillColor: loc.isHome ? '#ef4444' : '#3b82f6',
            color: '#fff',
            weight: 2,
            opacity: 1,
            fillOpacity: 0.8
          }).addTo(map);

          // Create popup content
          let popupContent = `
            <div class="map-popup">
              <h4>${loc.name}</h4>
              <p><strong>${loc.city}</strong></p>
              ${loc.collaborator ? `<p>Collaborator: ${loc.collaborator}</p>` : ''}
              <p class="text-sm text-gray-600 mt-2">Publications:</p>
              <ul class="text-xs text-gray-700 ml-4 list-disc">
                ${loc.publications.map(pub => `<li>${pub}</li>`).join('')}
              </ul>
              <span class="publication-count">${loc.publications.length} publication${loc.publications.length > 1 ? 's' : ''}</span>
            </div>
          `;

          marker.bindPopup(popupContent);

          // Draw lines from home locations to collaborators
          if (!loc.isHome) {
            const homeCoords = locations.find(l => l.isHome && l.name === 'Los Alamos National Laboratory (LANL)');
            if (homeCoords) {
              L.polyline([homeCoords.coords, loc.coords], {
                color: '#3b82f6',
                weight: 2,
                opacity: 0.4,
                dashArray: '5, 10'
              }).addTo(map);
            }
          }
        });

        // Add legend with reset button
        const legend = L.control({ position: 'bottomright' });
        legend.onAdd = function () {
          const div = L.DomUtil.create('div', 'map-legend');
          div.style.cssText = 'background: white; padding: 12px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.15);';
          div.innerHTML = `
            <div style="font-size: 14px; font-weight: 600; margin-bottom: 8px; color: #1e293b;">Legend</div>
            <div style="display: flex; align-items: center; margin-bottom: 4px;">
              <div style="width: 12px; height: 12px; border-radius: 50%; background: #ef4444; margin-right: 8px;"></div>
              <span style="font-size: 12px; color: #64748b;">Home Institution</span>
            </div>
            <div style="display: flex; align-items: center; margin-bottom: 12px;">
              <div style="width: 8px; height: 8px; border-radius: 50%; background: #3b82f6; margin-right: 8px;"></div>
              <span style="font-size: 12px; color: #64748b;">Collaborator</span>
            </div>
            <button onclick="resetCollaborationMap()" style="width: 100%; background: #3b82f6; color: white; font-size: 12px; padding: 6px 12px; border: none; border-radius: 4px; cursor: pointer; transition: background 0.3s;" onmouseover="this.style.background='#2563eb'" onmouseout="this.style.background='#3b82f6'">
              Reset View
            </button>
          `;
          return div;
        };
        legend.addTo(map);

        // Store map instance globally for reset button
        window.collaborationMapInstance = map;

      } catch (error) {
        console.error('Error initializing collaboration map:', error);
      }
    }

    /* ---------- Reset Collaboration Map ---------- */
    window.resetCollaborationMap = function () {
      if (window.collaborationMapInstance) {
        window.collaborationMapInstance.setView([39.8283, -98.5795], 4);
      }
    };

    /* ---------- Project Filtering ---------- */
    window.filterProjects = function (category) {
      const projects = document.querySelectorAll('.project-item');
      const buttons = document.querySelectorAll('.filter-btn');

      // Update active button
      buttons.forEach(btn => {
        btn.classList.toggle('active', btn.dataset.category === category);
      });

      // Filter projects
      projects.forEach(project => {
        if (category === 'all' || project.dataset.category === category) {
          project.classList.remove('filtered-out');
        } else {
          project.classList.add('filtered-out');
        }
      });
    };

    /* ---------- Update initAll to include new features ---------- */
    // Note: initCollaborationMap will be called separately after Leaflet loads
    window.initCollaborationMap = initCollaborationMap;

    /* ---------- Exports (if needed) ---------- */
    window.portfolioEnhancements = {
      initParticles,
      initTypedText,
      initAOS,
      initCounters,
      initProgressBar,
      initSmoothScroll,
      initCollaborationMap,
    };
  })();
}
