// Enhanced Scripts for Ajay Khanna's Portfolio
// Phase 1 Implementation: Hero Section, Particles, Typed.js, Counters, AOS

// Wait for DOM to be fully loaded
document.addEventListener('DOMContentLoaded', function() {

  // Initialize all components
  initParticles();
  initTypedText();
  initAOS();
  initCounters();
  initProgressBar();
  initSmoothScroll();

});

// 1. Initialize Particles.js for animated background
function initParticles() {
  if (typeof particlesJS === 'undefined') {
    console.log('Particles.js not loaded');
    return;
  }

  particlesJS('particles-js', {
    particles: {
      number: {
        value: 80,
        density: {
          enable: true,
          value_area: 800
        }
      },
      color: {
        value: ['#3b82f6', '#10b981', '#14b8a6', '#8b5cf6']
      },
      shape: {
        type: 'circle',
        stroke: {
          width: 0,
          color: '#000000'
        }
      },
      opacity: {
        value: 0.5,
        random: true,
        anim: {
          enable: true,
          speed: 1,
          opacity_min: 0.1,
          sync: false
        }
      },
      size: {
        value: 3,
        random: true,
        anim: {
          enable: true,
          speed: 2,
          size_min: 0.1,
          sync: false
        }
      },
      line_linked: {
        enable: true,
        distance: 150,
        color: '#ffffff',
        opacity: 0.2,
        width: 1
      },
      move: {
        enable: true,
        speed: 2,
        direction: 'none',
        random: false,
        straight: false,
        out_mode: 'out',
        bounce: false,
        attract: {
          enable: false,
          rotateX: 600,
          rotateY: 1200
        }
      }
    },
    interactivity: {
      detect_on: 'canvas',
      events: {
        onhover: {
          enable: true,
          mode: 'grab'
        },
        onclick: {
          enable: true,
          mode: 'push'
        },
        resize: true
      },
      modes: {
        grab: {
          distance: 140,
          line_linked: {
            opacity: 0.5
          }
        },
        push: {
          particles_nb: 4
        },
        remove: {
          particles_nb: 2
        }
      }
    },
    retina_detect: true
  });
}

// 2. Initialize Typed.js for dynamic typing animation
function initTypedText() {
  if (typeof Typed === 'undefined') {
    console.log('Typed.js not loaded');
    return;
  }

  const typedElement = document.getElementById('typed-output');
  if (typedElement) {
    new Typed('#typed-output', {
      strings: [
        'Drug Discovery Expert',
        'Energy Transfer Researcher',
        'Battery Materials Scientist',
        'Green Energy Innovator',
        'ML for Chemistry Enthusiast',
        'QM/MM Specialist'
      ],
      typeSpeed: 50,
      backSpeed: 30,
      backDelay: 2000,
      loop: true,
      showCursor: true,
      cursorChar: '|',
      autoInsertCss: true
    });
  }
}

// 3. Initialize AOS (Animate On Scroll)
function initAOS() {
  if (typeof AOS === 'undefined') {
    console.log('AOS not loaded');
    return;
  }

  AOS.init({
    duration: 1000,
    once: false,
    mirror: true,
    offset: 100,
    easing: 'ease-in-out',
    anchorPlacement: 'top-bottom'
  });
}

// 4. Animated Counter for Statistics
function initCounters() {
  const counters = document.querySelectorAll('.counter');

  if (counters.length === 0) return;

  // Counter animation function
  function animateCounter(counter) {
    const target = parseInt(counter.getAttribute('data-target'));
    const duration = 2000; // 2 seconds
    const steps = 60;
    const increment = target / steps;
    let current = 0;

    const timer = setInterval(() => {
      current += increment;
      if (current >= target) {
        counter.textContent = target + '+';
        clearInterval(timer);
      } else {
        counter.textContent = Math.floor(current);
      }
    }, duration / steps);
  }

  // Intersection Observer to trigger counter when visible
  const observerOptions = {
    threshold: 0.5,
    rootMargin: '0px'
  };

  const observer = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      if (entry.isIntersecting && !entry.target.classList.contains('counted')) {
        entry.target.classList.add('counted');
        animateCounter(entry.target);
      }
    });
  }, observerOptions);

  counters.forEach(counter => {
    observer.observe(counter);
  });
}

// 5. Reading Progress Bar
function initProgressBar() {
  const progressBar = document.getElementById('progress-bar');

  if (!progressBar) return;

  window.addEventListener('scroll', () => {
    const windowHeight = window.innerHeight;
    const documentHeight = document.documentElement.scrollHeight;
    const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
    const scrollPercent = (scrollTop / (documentHeight - windowHeight)) * 100;

    progressBar.style.width = Math.min(scrollPercent, 100) + '%';
  });
}

// 6. Smooth Scroll for Anchor Links
function initSmoothScroll() {
  document.querySelectorAll('a[href^="#"]').forEach(anchor => {
    anchor.addEventListener('click', function (e) {
      const href = this.getAttribute('href');

      // Don't prevent default for just "#"
      if (href === '#') return;

      e.preventDefault();

      const target = document.querySelector(href);
      if (target) {
        const navHeight = document.querySelector('nav').offsetHeight;
        const targetPosition = target.getBoundingClientRect().top + window.pageYOffset - navHeight;

        window.scrollTo({
          top: targetPosition,
          behavior: 'smooth'
        });
      }
    });
  });
}

// 7. Update Navigation Toggle (Keep existing functionality)
const navToggle = document.getElementById('nav-toggle');
const navContent = document.getElementById('nav-content');

if (navToggle && navContent) {
  navToggle.addEventListener('click', function () {
    navContent.classList.toggle('hidden');
  });
}

// 8. Project Slider Functionality (Keep existing)
let currentPosition = 0;
const slider = document.getElementById('projectSlider');

if (slider) {
  const slideWidth = slider.children[0]?.offsetWidth || 0;
  const totalSlides = slider.children.length;

  window.slide = function(direction) {
    currentPosition = Math.max(Math.min(currentPosition + direction, 0), -totalSlides + 1);
    slider.style.transform = `translateX(${currentPosition * slideWidth}px)`;
  };
}

// 9. Tutorial Slider Functionality (Keep existing)
let currentPositionTutorials = 0;
const tutorialSlider = document.getElementById('tutorialSlider');

if (tutorialSlider) {
  const tutorialSlideWidth = tutorialSlider.children[0]?.offsetWidth || 0;
  const totalTutorialSlides = tutorialSlider.children.length;

  window.slideTutorials = function(direction) {
    currentPositionTutorials = Math.max(Math.min(currentPositionTutorials + direction, 0), -totalTutorialSlides + 1);
    tutorialSlider.style.transform = `translateX(${currentPositionTutorials * tutorialSlideWidth}px)`;
  };
}

// 10. Handle Window Resize for Sliders
window.addEventListener('resize', function() {
  // Update slider positions on resize
  if (slider && slider.children.length > 0) {
    const newSlideWidth = slider.children[0].offsetWidth;
    slider.style.transform = `translateX(${currentPosition * newSlideWidth}px)`;
  }

  if (tutorialSlider && tutorialSlider.children.length > 0) {
    const newTutorialWidth = tutorialSlider.children[0].offsetWidth;
    tutorialSlider.style.transform = `translateX(${currentPositionTutorials * newTutorialWidth}px)`;
  }
});

// 11. Performance Optimization: Lazy Loading Images
function initLazyLoading() {
  if ('IntersectionObserver' in window) {
    const imageObserver = new IntersectionObserver((entries, observer) => {
      entries.forEach(entry => {
        if (entry.isIntersecting) {
          const img = entry.target;
          img.src = img.dataset.src;
          img.classList.remove('lazy');
          observer.unobserve(img);
        }
      });
    });

    document.querySelectorAll('img.lazy').forEach(img => {
      imageObserver.observe(img);
    });
  }
}

// Initialize lazy loading after DOM is ready
initLazyLoading();

// 12. Add active state to navigation items based on scroll position
function updateActiveNav() {
  const sections = document.querySelectorAll('section[id]');
  const navLinks = document.querySelectorAll('nav a[href^="#"]');

  let current = '';

  sections.forEach(section => {
    const sectionTop = section.offsetTop;
    const sectionHeight = section.clientHeight;
    if (pageYOffset >= (sectionTop - 200)) {
      current = section.getAttribute('id');
    }
  });

  navLinks.forEach(link => {
    link.classList.remove('active');
    if (link.getAttribute('href') === '#' + current) {
      link.classList.add('active');
    }
  });
}

window.addEventListener('scroll', updateActiveNav);

// 13. Console Message (Fun Easter Egg)
console.log('%c👋 Hello! ', 'font-size: 20px; font-weight: bold; color: #3b82f6;');
console.log('%cInterested in my work? Let\'s collaborate!', 'font-size: 14px; color: #10b981;');
console.log('%c🔬 Computational Chemistry | 💊 Drug Discovery | ⚡ Energy Solutions', 'font-size: 12px; color: #8b5cf6;');

// Export functions for potential use in other scripts
window.portfolioEnhancements = {
  initParticles,
  initTypedText,
  initAOS,
  initCounters,
  initProgressBar,
  initSmoothScroll
};
