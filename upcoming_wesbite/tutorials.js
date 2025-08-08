/**
 * @file Manages the tutorials page functionality.
 * @summary This script fetches tutorial data from a manifest file,
 * dynamically renders them on the page, and provides search and filtering
 * capabilities. It also includes a modal to view the full tutorial content.
 *
 * @description
 * The script performs the following actions upon page load:
 * 1. Fetches a `manifest.json` file which lists all available tutorial HTML files.
 * 2. For each file in the manifest, it fetches the HTML content.
 * 3. It parses each HTML file to extract metadata (title, description, image, tags)
 *    from meta tags and the main content from the `<body>`.
 * 4. It populates a grid with "tutorial cards" and displays the first tutorial
 *    as a "featured" item.
 * 5. It sets up event listeners for a search bar to filter tutorials in real-time
 *    based on title, description, or tags.
 * 6. It handles clicks on tutorial cards to open a modal window displaying the
 *    full tutorial content.
 * 7. It includes functionality to close the modal.
 * 8. It handles errors during data fetching and displays a "no results" message
 *    if tutorials cannot be loaded or if a search yields no results.
 */

document.addEventListener("DOMContentLoaded", () => {
  // Element selectors
  const tutorialsGrid = document.getElementById("tutorials-grid");
  const featuredTutorial = document.getElementById("featured-tutorial");
  const searchBar = document.getElementById("search-bar");
  const searchButton = document.getElementById("search-button");
  const noResults = document.getElementById("no-results");
  const modal = document.getElementById("tutorial-modal");
  const modalBody = document.getElementById("modal-body");
  const closeButton = document.querySelector(".close-button");

  let allTutorials = [];

  // --- Data Fetching and Initialization ---

  /**
   * Fetches the manifest file, then fetches and parses all listed tutorials.
   */
  fetch("tutorials/manifest.json")
    .then((response) => {
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
      return response.json();
    })
    .then((manifest) => {
      if (!manifest.tutorials || manifest.tutorials.length === 0) {
        displayTutorials([]); // Show no results if manifest is empty
        return;
      }
      const tutorialPromises = manifest.tutorials.map((tutorialFile) =>
        fetch(`./tutorials/${tutorialFile}`).then((res) => {
          if (!res.ok) {
            console.error(`Failed to load tutorial: ${tutorialFile}`);
            return null; // Return null for failed fetches to handle them later
          }
          return res.text();
        })
      );

      Promise.all(tutorialPromises).then((tutorialsData) => {
        // Filter out any null results from failed fetches and parse the successful ones
        allTutorials = tutorialsData
          .filter((data) => data !== null)
          .map((data) => {
            const parser = new DOMParser();
            const doc = parser.parseFromString(data, "text/html");

            // Safely extract metadata with fallbacks
            const title =
              doc
                .querySelector('meta[name="title"]')
                ?.getAttribute("content") || "No Title";
            const description =
              doc
                .querySelector('meta[name="description"]')
                ?.getAttribute("content") || "No description available.";
            const image =
              doc
                .querySelector('meta[name="image"]')
                ?.getAttribute("content") ||
              "https://placehold.co/600x400/E2E8F0/4A5568?text=No+Image";
            const tags =
              doc.querySelector('meta[name="tags"]')?.getAttribute("content") ||
              "";
            const body = doc.querySelector("body").innerHTML;

            return { title, description, image, tags, body };
          });

        displayTutorials(allTutorials);
        if (allTutorials.length > 0) {
          displayFeaturedTutorial(allTutorials[0]);
        }
      });
    })
    .catch((error) => {
      console.error("Failed to load tutorials manifest:", error);
      displayTutorials([]); // Show the "no results" message on error
    });

  // --- UI Rendering Functions ---

  /**
   * Renders the tutorial cards in the grid or shows the "no results" message.
   * @param {Array} tutorials - The array of tutorial objects to display.
   */
  function displayTutorials(tutorials) {
    tutorialsGrid.innerHTML = "";
    if (tutorials.length === 0) {
      tutorialsGrid.style.display = "none";
      noResults.style.display = "block";
      featuredTutorial.style.display = "none"; // Also hide featured tutorial
    } else {
      tutorialsGrid.style.display = "grid";
      noResults.style.display = "none";
      tutorials.forEach((tutorial) => {
        const card = document.createElement("div");
        card.className =
          "tutorial-card bg-white rounded-lg shadow-lg p-4 cursor-pointer transition transform hover:scale-105";
        card.innerHTML = `
                    <img src="${tutorial.image}" alt="${tutorial.title}" class="w-full h-48 object-cover mb-4 rounded-md" onerror="this.onerror=null;this.src='https://placehold.co/600x400/E2E8F0/4A5568?text=Image+Error';">
                    <h3 class="text-xl font-semibold mb-2">${tutorial.title}</h3>
                    <p class="text-gray-700 mb-4">${tutorial.description}</p>
                    <div class="tags text-sm text-gray-500"><strong>Tags:</strong> ${tutorial.tags}</div>
                `;
        card.addEventListener("click", () => openModal(tutorial));
        tutorialsGrid.appendChild(card);
      });
    }
  }

  /**
   * Renders the featured tutorial section.
   * @param {Object} tutorial - The tutorial object to feature.
   */
  function displayFeaturedTutorial(tutorial) {
    if (!tutorial) return;
    featuredTutorial.innerHTML = `
            <h2 class="text-3xl font-bold mb-4">Featured Tutorial</h2>
            <div class="featured-card flex flex-col md:flex-row bg-white rounded-lg shadow-lg overflow-hidden">
                <img src="${tutorial.image}" alt="${
      tutorial.title
    }" class="w-full md:w-1/3 object-cover" onerror="this.onerror=null;this.src='https://placehold.co/600x400/E2E8F0/4A5568?text=Image+Error';">
                <div class="w-full md:w-2/3 p-6">
                    <h3 class="text-2xl font-bold mb-2">${tutorial.title}</h3>
                    <p class="text-gray-700 mb-4">${tutorial.description}</p>
                    <button class="bg-blue-500 text-white px-6 py-2 rounded hover:bg-blue-600 transition-colors" onclick="openModalWithData('${btoa(
                      JSON.stringify(tutorial)
                    )}')">Read More</button>
                </div>
            </div>
        `;
    featuredTutorial.style.display = "block";
  }

  // --- Event Handlers ---

  /**
   * Filters and displays tutorials based on the search bar's value.
   */
  function handleSearch() {
    const searchTerm = searchBar.value.toLowerCase();
    const filteredTutorials = allTutorials.filter(
      (tutorial) =>
        tutorial.title.toLowerCase().includes(searchTerm) ||
        tutorial.description.toLowerCase().includes(searchTerm) ||
        tutorial.tags.toLowerCase().includes(searchTerm)
    );
    displayTutorials(filteredTutorials);
    if (filteredTutorials.length > 0) {
      displayFeaturedTutorial(filteredTutorials[0]);
    }
  }

  searchButton.addEventListener("click", handleSearch);
  searchBar.addEventListener("keyup", (e) => {
    if (e.key === "Enter") {
      handleSearch();
    }
  });

  /**
   * Opens the modal with the content of the selected tutorial.
   * @param {Object} tutorial - The tutorial object to display in the modal.
   */
  function openModal(tutorial) {
    modalBody.innerHTML = tutorial.body;
    modal.style.display = "block";
    // Re-run highlight.js on the new content in the modal
    document.querySelectorAll("#modal-body pre code").forEach((block) => {
      hljs.highlightElement(block);
    });
  }

  // A globally accessible version of openModal to be called from the featured button's onclick attribute.
  // It decodes the base64 data string back into a tutorial object.
  window.openModalWithData = (base64Data) => {
    try {
      const tutorial = JSON.parse(atob(base64Data));
      openModal(tutorial);
    } catch (e) {
      console.error("Failed to parse tutorial data for modal:", e);
    }
  };

  // Close modal events
  closeButton.addEventListener("click", () => {
    modal.style.display = "none";
  });

  window.addEventListener("click", (e) => {
    if (e.target == modal) {
      modal.style.display = "none";
    }
  });
});
