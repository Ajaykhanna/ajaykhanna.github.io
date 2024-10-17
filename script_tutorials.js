document.addEventListener('DOMContentLoaded', function () {
    hljs.highlightAll();

    const tutorialCards = document.getElementsByClassName('tutorial-card');
    const tutorialSection = document.getElementById('tutorial-section');

    for (let i = 0; i < tutorialCards.length; i++) {
        const card = tutorialCards[i];
        const link = card.querySelector('a');
        link.addEventListener('click', function (event) {
            event.preventDefault();
            const tutorialTitle = card.querySelector('h3').innerText;
            loadTutorial(tutorialTitle);
        });
    }

    function loadTutorial(title) {
        // Clear any existing tutorial content
        while (tutorialSection.firstChild) {
            tutorialSection.removeChild(tutorialSection.firstChild);
        }

        let detailedContent = '';
        switch (title) {
            case 'Quantum Mechanics Basics':
                detailedContent = `
                    <p style="width: 700px; display: inline-block;">In this tutorial, you will learn about the basics of quantum mechanics, including wave functions, the Schrödinger equation, and operators. We will also explore practical implementations using Gaussian software to solve real-world quantum problems.</p>
                    <img src="quantum_example.png" alt="Quantum Mechanics Example" class="w-auto mb-4">
                    <p style="width: 700px; display: inline-block;">Example topics covered include:</p>
                    <ul>
                        <li>Introduction to Quantum States</li>
                        <li>Wave-Particle Duality</li>
                        <li>Setting up Gaussian for simple quantum systems</li>
                    </ul>
                    <p style="width: 500px; display: inline-block;">Here is an example of a Python script used in quantum simulations:</p>
                    <div class="code-box">
                        <pre><code class="language-python">
# Importing necessary modules
import numpy as np

# Define the Hamiltonian matrix
H = np.array([[1, 0.1], [0.1, 2]])

# Calculate eigenvalues and eigenvectors
eigenvalues, eigenvectors = np.linalg.eigh(H)

# Print the results
print("Eigenvalues:", eigenvalues)
print("Eigenvectors:", eigenvectors)
                        </code></pre>
                    </div>`;
                break;

            case 'Molecular Dynamics with GROMACS':
                detailedContent = `
                    <p style="width: 700px; display: inline-block;">This tutorial provides a comprehensive guide to running molecular dynamics (MD) simulations using GROMACS. We will cover all aspects, from installation and configuration to analyzing simulation trajectories.</p>
                    <img src="gromacs_example.png" alt="Molecular Dynamics Example" class="w-auto mb-4">
                    <p style="width: 700px; display: inline-block;">Key topics include:</p>
                    <ul>
                        <li>Setting up a simulation box</li>
                        <li>Energy minimization and equilibration</li>
                        <li>Production run and data analysis</li>
                    </ul>
                    <p style="width: 500px; display: inline-block;">Below is an example of a Python script for analyzing GROMACS output:</p>
                    <div class="code-box">
                        <pre><code class="language-python">
# Import MDAnalysis to process GROMACS trajectory
import MDAnalysis as mda

# Load the trajectory and topology files
u = mda.Universe('topol.tpr', 'traj.xtc')

# Select all atoms
atoms = u.select_atoms('all')

# Iterate over each frame and calculate the center of mass
for ts in u.trajectory:
    print("Center of Mass:", atoms.center_of_mass())
                        </code></pre>
                    </div>`;
                break;

            case 'Hybrid QM/MM Simulations':
                detailedContent = `
                    <p style="width: 700px; display: inline-block;">Hybrid QM/MM simulations combine quantum mechanics and molecular mechanics to study large molecular systems. This tutorial will help you understand how to implement QM/MM using software like TeraChem and Amber.</p>
                    <img src="qm_mm_example.png" alt="QM/MM Example" class="w-auto mb-4">
                    <p style="width: 700px; display: inline-block;">Topics covered include:</p>
                    <ul>
                        <li>Choosing the QM and MM regions</li>
                        <li>Energy calculations using hybrid methods</li>
                        <li>Practical tips for efficient QM/MM simulations</li>
                    </ul>
                    <p style="width: 500px; display: inline-block;">Here is an example Python code for partitioning QM and MM regions:</p>
                    <div class="code-box">
                        <pre><code class="language-python">
# Define QM and MM regions
qm_atoms = [1, 2, 3]  # List of QM atom indices
mm_atoms = [4, 5, 6]  # List of MM atom indices

# Calculate QM and MM interactions
def calculate_qm_energy(qm_atoms):
    # Placeholder function for QM calculations
    return sum(qm_atoms) * 0.1

def calculate_mm_energy(mm_atoms):
    # Placeholder function for MM calculations
    return sum(mm_atoms) * 0.05

# Print the energies
print("QM Energy:", calculate_qm_energy(qm_atoms))
print("MM Energy:", calculate_mm_energy(mm_atoms))
                        </code></pre>
                    </div>`;
                break;

            default:
                detailedContent = `
                    <p style="width: 700px; display: inline-block;">Here is a more detailed explanation of the topic, including examples, figures, and additional resources to help you understand the concepts better. You can also find practical applications and step-by-step instructions to follow along with the tutorial.</p>`;
        }

        const tutorialContainer = document.createElement('div');
        tutorialContainer.className = 'tutorial-content';
        tutorialContainer.innerHTML = `
            <div class="close-btn">&times; Close</div>
            <button class="maximize-btn">⤢ Maximize</button>
            <h2 class="text-2xl font-bold mb-4">${title}</h2>
            ${detailedContent}
        `;
        tutorialSection.appendChild(tutorialContainer);
        tutorialSection.style.display = 'flex';
        tutorialSection.style.overflowX = 'auto';
        tutorialSection.style.whiteSpace = 'wrap';

        // Close button functionality
        tutorialContainer.querySelector('.close-btn').addEventListener('click', function () {
            tutorialSection.removeChild(tutorialContainer);
            if (tutorialSection.children.length === 0) {
                tutorialSection.style.display = 'none';
            }
        });

        // Maximize button functionality
        tutorialContainer.querySelector('.maximize-btn').addEventListener('click', function () {
            tutorialContainer.classList.toggle('maximized');
            if (tutorialContainer.classList.contains('maximized')) {
                tutorialContainer.style.position = 'fixed';
                tutorialContainer.style.top = '0';
                tutorialContainer.style.left = '0';
                tutorialContainer.style.width = '100%';
                tutorialContainer.style.height = '100%';
                tutorialContainer.style.zIndex = '1000';
                tutorialContainer.style.backgroundColor = 'white';
                tutorialContainer.style.overflow = 'auto';
            } else {
                tutorialContainer.style.position = 'static';
                tutorialContainer.style.width = '600px';
                tutorialContainer.style.height = '400px';
                tutorialContainer.style.zIndex = '1';
                tutorialContainer.style.backgroundColor = '#f8fafc';
                tutorialContainer.style.overflow = 'hidden';
            }
        });
    }
});
