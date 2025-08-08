// scripts.js

/**
 * Toggle navigation bar visibility for mobile devices.
 */
document.getElementById('nav-toggle').addEventListener('click', function () {
    const navContent = document.getElementById('nav-content');
    navContent.classList.toggle('hidden');
});

/**
 * Slides project cards left or right.
 * @param {number} direction - The direction to slide. Positive for right, negative for left.
 */
let currentPosition = 0;
const slider = document.getElementById('projectSlider');
const slideWidth = slider ? slider.children[0].offsetWidth : 0; // Ensure slider exists
const totalSlides = slider ? slider.children.length : 0;

function slide(direction) {
    if (slider) {
        currentPosition = Math.max(Math.min(currentPosition + direction, 0), -totalSlides + 1);
        slider.style.transform = `translateX(${currentPosition * slideWidth}px)`;
    }
}

/**
 * Slides tutorial cards left or right.
 * @param {number} direction - The direction to slide. Positive for right, negative for left.
 */
let currentPositionTutorials = 0;
const tutorialSlider = document.getElementById('tutorialSlider');
const tutorialSlideWidth = tutorialSlider ? tutorialSlider.children[0].offsetWidth : 0; // Ensure slider exists
const totalTutorialSlides = tutorialSlider ? tutorialSlider.children.length : 0;

function slideTutorials(direction) {
    if (tutorialSlider) {
        currentPositionTutorials = Math.max(Math.min(currentPositionTutorials + direction, 0), -totalTutorialSlides + 1);
        tutorialSlider.style.transform = `translateX(${currentPositionTutorials * tutorialSlideWidth}px)`;
    }
}
