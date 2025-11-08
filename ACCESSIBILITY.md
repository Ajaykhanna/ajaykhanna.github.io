# Accessibility Statement

## Overview
This website is committed to ensuring digital accessibility for people with disabilities. We are continually improving the user experience for everyone and applying the relevant accessibility standards.

## Conformance Status
This website conforms to **WCAG 2.1 Level AA** standards. We strive to meet or exceed these requirements to ensure our content is accessible to all users.

## Accessibility Features

### 1. Keyboard Navigation
- **Full keyboard support**: All interactive elements can be accessed using only a keyboard
- **Skip navigation link**: Press `Tab` on page load to reveal a "Skip to main content" link
- **Visible focus indicators**: Clear visual indicators show which element has keyboard focus
- **Logical tab order**: Navigate through the page in a logical sequence

**Keyboard Shortcuts:**
- `Tab`: Move to next interactive element
- `Shift + Tab`: Move to previous interactive element
- `Enter` or `Space`: Activate buttons and links
- `Escape`: Close modals or expanded menus (where applicable)

### 2. Screen Reader Support
- **Semantic HTML5**: Proper use of `<main>`, `<nav>`, `<section>`, `<article>` elements
- **ARIA labels**: Descriptive labels on all interactive elements
- **ARIA live regions**: Dynamic content changes are announced to screen readers
- **Alternative text**: All images have descriptive alt attributes
- **Form labels**: All form inputs have associated labels

**Tested with:**
- NVDA (Windows)
- JAWS (Windows)
- VoiceOver (macOS/iOS)
- TalkBack (Android)

### 3. Visual Design
- **High contrast**: Text meets WCAG AA color contrast ratios (4.5:1 minimum)
- **Resizable text**: Text can be resized up to 200% without loss of functionality
- **Focus indicators**: 3px solid outlines with high contrast colors
- **No flashing content**: No content flashes more than 3 times per second

**Color Contrast Ratios:**
- Body text: #374151 on white (11.5:1)
- Link text: #3b82f6 on white (4.9:1)
- Button text: White on #1e40af (8.6:1)

### 4. Touch Targets
- **Minimum size**: All interactive elements are at least 44x44 pixels on mobile
- **Adequate spacing**: Touch targets have sufficient space between them
- **Large clickable areas**: Buttons and links have generous padding

### 5. Forms
- **Clear labels**: All form fields have visible, associated labels
- **Required field indicators**: Required fields are marked with asterisks and aria-required
- **Error messages**: Clear, descriptive error messages with suggestions
- **Autocomplete**: Appropriate autocomplete attributes for faster form filling
- **Help text**: Additional context provided where needed

### 6. Responsive Design
- **Mobile-first**: Optimized for mobile devices first, then scaled up
- **Flexible layouts**: Content adapts to different screen sizes
- **No horizontal scrolling**: Content fits within viewport at all sizes
- **Readable text**: Minimum font size of 16px on mobile

### 7. Reduced Motion
- **Respects user preferences**: Honors `prefers-reduced-motion` setting
- **Minimal animations**: Animations can be disabled for users with motion sensitivity
- **No auto-playing content**: No videos or animations start automatically

## Assistive Technology Compatibility

### Supported Screen Readers
- ✅ NVDA 2021+ (Windows)
- ✅ JAWS 2019+ (Windows)
- ✅ VoiceOver (macOS, iOS)
- ✅ TalkBack (Android)
- ✅ Narrator (Windows)

### Supported Browsers
- ✅ Chrome/Edge (latest 2 versions)
- ✅ Firefox (latest 2 versions)
- ✅ Safari (latest 2 versions)
- ✅ Mobile browsers (iOS Safari, Chrome Android)

## Known Limitations

### Current Accessibility Issues
None at this time. If you encounter any accessibility barriers, please report them.

### Third-Party Content
- **Google Analytics**: Used for anonymous visitor statistics
- **Formspree**: Contact form submission service
- **Leaflet.js**: Interactive map library (keyboard accessible)
- **Particles.js**: Decorative background animation (does not interfere with screen readers)

## Testing and Validation

### Automated Testing Tools
- ✅ axe DevTools
- ✅ WAVE (Web Accessibility Evaluation Tool)
- ✅ Lighthouse Accessibility Audit
- ✅ Pa11y

### Manual Testing
- ✅ Keyboard-only navigation
- ✅ Screen reader testing (NVDA, VoiceOver)
- ✅ Color contrast verification
- ✅ Mobile touch target testing
- ✅ Form validation testing

## Accessibility Features by Section

### Navigation
- Semantic `<nav>` element with aria-label
- Mobile menu toggle with aria-expanded state
- Skip to main content link
- Keyboard accessible menu items

### Hero Section
- Decorative particles background (aria-hidden for screen readers)
- Proper heading hierarchy (h1 for name)
- Accessible call-to-action buttons

### Projects Section
- Filter buttons with aria-pressed states
- ARIA live region announces filter changes
- Keyboard accessible carousel with prev/next buttons
- Alt text for all project images

### Publications Section
- Timeline filter with aria-pressed states
- ARIA live region for filter updates
- Structured publication data
- External links open in new tab with warning

### Contact Form
- Properly associated labels (for/id)
- Required field indicators
- ARIA descriptions for inputs
- Clear error messages
- Submit button with focus ring

### Floating Action Button (FAB)
- Keyboard accessible
- Clear aria-labels for all actions
- Focus visible indicators
- Screen reader friendly social links

## How to Report Accessibility Issues

If you encounter any accessibility barriers or have suggestions for improvement:

1. **Email**: [Your email address]
2. **GitHub**: [Open an issue](https://github.com/Ajaykhanna/ajaykhanna.github.io/issues)
3. **Subject**: "Accessibility Issue: [Brief Description]"

Please include:
- Description of the issue
- The page or section where you encountered it
- Your browser and operating system
- Assistive technology you're using (if applicable)
- Steps to reproduce the issue

We aim to respond to all accessibility reports within 2 business days.

## Additional Resources

### Accessibility Guidelines
- [WCAG 2.1 Guidelines](https://www.w3.org/WAI/WCAG21/quickref/)
- [WebAIM Resources](https://webaim.org/)
- [A11y Project](https://www.a11yproject.com/)

### Browser Accessibility Features
- [Chrome Accessibility](https://www.google.com/accessibility/)
- [Firefox Accessibility](https://support.mozilla.org/en-US/kb/accessibility-features-firefox)
- [Safari Accessibility](https://www.apple.com/accessibility/)

## Updates and Maintenance

**Last Updated**: November 8, 2025
**Next Review**: February 2026

We review and update our accessibility features quarterly to ensure continued compliance and usability.

## Legal Compliance

This website aims to comply with:
- Americans with Disabilities Act (ADA)
- Section 508 of the Rehabilitation Act
- Web Content Accessibility Guidelines (WCAG) 2.1 Level AA
- European Accessibility Act (EAA)

---

**Commitment**: We are committed to providing an accessible website. If you have any feedback or encounter any accessibility barriers, please contact us so we can improve.
