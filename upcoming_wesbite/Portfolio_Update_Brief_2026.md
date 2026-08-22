# Portfolio Website Update Brief — ajaykhanna.github.io

**For:** a fresh Claude Code session (no prior context). This document is self-contained.
**Goal:** review-driven update of Ajay Khanna's personal portfolio site. Make the verified fixes,
reposition the content around his strongest provable expertise, and apply polish — without breaking
the existing look and feel.

> Owner preferences (must follow):
> 1. **Contact = a simple Google Form.** Replace the existing (broken) form with a link/button to a Google Form. Ajay will give you the Google Form URL — ask for it; until then use the placeholder `{{GOOGLE_FORM_URL}}`.
> 2. **No public CV/resume download.** Privacy + anti-spam: people must reach out via the contact Google Form first, then Ajay shares the CV manually. Remove every public resume link/file and replace it with a "CV available on request" call-to-action that points to the Google Form.

---

## 0. Setup & ground facts about the site

- **Live URL:** https://ajaykhanna.github.io/
- **Repo:** `Ajaykhanna/ajaykhanna.github.io` (GitHub Pages, user site, served from the default branch root).
  - Clone: `gh repo clone Ajaykhanna/ajaykhanna.github.io` (or `git clone https://github.com/Ajaykhanna/ajaykhanna.github.io`).
- **Stack (hand-coded, no static-site generator):**
  - Single page: **`index.html`** (~115 KB, all sections inline).
  - **Tailwind CSS via the Play CDN** (`cdn.tailwindcss.com`) — dev CDN, not a build.
  - `particles.js` (hero background) and `typed.js` (animated headline) via jsDelivr.
  - A separate **`resume.html`** page exists (currently linked publicly — see Task B).
- **Last updated:** 2026-03-04 (commit "Updated Publications") — i.e. stale by several months.
- **Local preview:** just open `index.html` in a browser, or `python -m http.server` in the repo and visit `localhost:8000`.
- **Nav/sections present (in current order):** Hero (stats) → Research Impact Areas (4 domains) → About → Projects → Technical Expertise → Experience → Featured Publications → Global Collaboration → Timeline → Certifications → Tutorials → Resume → Contact.
- **Work style:** before editing, read `index.html` top-to-bottom so you know the existing Tailwind class patterns and section anchors (`#about`, `#projects`, `#skills`, `#experience`, `#tutorials`, `#resume`, `#contact`). Match the existing visual language; do not introduce a new framework.

---

## 1. IMMEDIATE FIXES (bugs / correctness) — do these first

### Task A — Fix the broken contact form (→ Google Form)
- **Problem (verified):** the contact `<form>` posts to `https://formspree.io/f/your_form_id` — a literal placeholder. **Every submission silently fails.**
- **Do:** remove the `<form>...</form>` and its name/email/message inputs. Replace with a clean Contact card that contains a prominent button:
  ```html
  <a href="{{GOOGLE_FORM_URL}}" target="_blank" rel="noopener"
     class="<match existing primary-button Tailwind classes>">Get in touch</a>
  ```
  Keep the existing social/profile links (GitHub, LinkedIn, Google Scholar, ORCID, ResearchGate, X) next to it. Optionally embed the form with an `<iframe>` instead of a button — but a button link is simpler and is what the owner asked for; prefer the button.
- **Ask the owner** for `{{GOOGLE_FORM_URL}}` and substitute it. Do not invent one.
- Update the Contact copy to set expectations, e.g. "The fastest way to reach me is this short form."

### Task B — Remove the public resume; gate the CV behind the form
- **Problem:** the site links a public resume (`#resume` section → `./resume.html`), and recruiters/visitors can reach it directly. Owner wants the CV **shared manually on request only**.
- **Do:**
  1. In the nav and the Resume section, **remove the public resume link/download.** Repurpose the section (keep the `#resume` anchor or rename to `#cv`) into a short call-to-action: e.g. "**CV available on request** — please reach out via the contact form and I'll be glad to share it." Button → `{{GOOGLE_FORM_URL}}`.
  2. **Make `resume.html` non-public:** the cleanest option is to **delete `resume.html`** from the repo (since the CV will be shared privately). If the owner wants to keep the file for his own use, instead add `<meta name="robots" content="noindex,nofollow">` to it AND remove every link to it — but note a GitHub Pages file is still reachable if the URL is guessed, so **deletion is the privacy-safe choice**. Confirm with the owner which he prefers; default to deletion.
  3. Ensure no `resume.pdf`/CV PDF is added anywhere (there is none today — keep it that way).
- **Acceptance:** no path on the site downloads or renders a CV; the only route to the CV is the contact form.

### Task C — Fix the contradictory citation counts
- **Problem (verified):** the Nature Communications paper is shown as **"16+ citations"** in the Experience section and **"22 citations"** in the Featured Publications section. Other hardcoded counts ("50 Citations" hero stat, per-paper counts) will also drift.
- **Do (pick one approach, prefer the durable one):**
  - **Preferred:** stop hardcoding citation counts. Replace the hero "Citations" stat and per-paper counts with a **live Google Scholar badge/link** (`https://scholar.google.com/citations?user=qJM0sOIAAAAJ`) and a single "Citations: see Google Scholar" link, so numbers never go stale or disagree.
  - **Minimum:** reconcile to one correct number everywhere (use the current Scholar value at edit time) and remove the duplicate/conflicting mention.

### Task D — De-stale the hero stats and dates
- The hero shows hardcoded "7 Publications, 50 Citations, 4 Open Source Tools, 2000 Articles Downloaded."
  - "**2000 Articles Downloaded**" is a vague, unsourced metric — remove it or replace with a concrete, verifiable signal (e.g. "Nature Communications author", "Journal cover", "MLCM-26 co-organizer").
  - Update publication/tool counts to current and keep them consistent with the Publications and Projects sections.
- **Mark in-press work:** the "Covalent Control of Excitonic Interactions… (Nano Letters, 2026)" entry is future-dated; label it "**in press / accepted, 2026**" so it does not look like a typo.

### Task E — (minor) Tailwind Play CDN
- The site loads `cdn.tailwindcss.com` (the Play CDN), which prints a "not for production" console warning and recompiles CSS in-browser (flash of unstyled content + slower first paint).
- **Optional improvement:** add a tiny build step to ship a static, purged CSS file (Tailwind CLI: `npx tailwindcss -i input.css -o dist.css --minify`) and link that instead. Low priority; only do it if time allows and it does not risk the layout.

---

## 2. STRUCTURAL / LAYOUT recommendations

**Keep the architecture** — a single-page Tailwind scroller with anchor nav is modern and appropriate. Do **not** rebuild it. Two targeted changes:

### Task F — Reorder to lead with proof
- Current order puts the aspirational "4 Research Impact Areas" grid **before** the publications. Flip the emphasis so credibility comes first:
  - **Suggested order:** Hero (sharper value prop — see Task H) → **Selected Publications + Open-Source Tools** (proof) → Research themes (2–3 honest pillars, Task H) → Experience → Timeline → Tutorials → Contact/CV-request.
- Keep all existing anchors working; update the nav order to match.

### Task G — (optional) Split the very long single page
- As content grows, consider moving Publications and Projects to their own pages (`/publications`, `/projects`) with a concise landing page. Benefits: scannability, SEO, deep-linking. Not required — only if the owner wants it. If you do, keep the same Tailwind styling and a shared header/nav.

---

## 3. CONTENT rewrite (highest-value work)

### Task H — Reposition around the strongest, *provable* expertise
- **Problem:** the site currently leads with four equally-weighted "cardinal domains": **Drug Discovery, Energy Transfer, Energy Storage (batteries/electrodes/Li-ion/solid-state), Green Energy (CO₂ reduction, water splitting, OPV).** The **Energy Storage** and **Green Energy** areas are aspirational — there are **no publications** backing them. Presenting them as established expertise is a credibility risk, and it buries his real differentiator.
- **His actual, publication-backed strengths** (lead with these): **machine-learned interatomic potentials (MLIPs) and ML-for-science** (GNN/message-passing models, active learning, uncertainty quantification, HPC), **excited-state / nonadiabatic molecular dynamics and spectroscopy**, and **computational drug discovery** (QM/MM free energy, structure-based design).
- **Do:** collapse the four tiles into **2–3 honest pillars**, each anchored to real publications/tools:
  1. **ML Interatomic Potentials & ML-for-Science** — HIP-NN / hippynn, GNN / message-passing networks, active learning, uncertainty quantification, GPU/HPC; tool: `mlip_benchmark`.
  2. **Excited-State & Nonadiabatic Dynamics / Spectroscopy** — his published line of work (Nano Lett, J. Phys. Chem. Lett. cover, J. Chem. Phys., Nature Communications).
  3. **Computational Drug Discovery** — QM/MM free energy, SBDD/CADD (Frontier Medicines internship: docking, MD, binding free energy).
- **Batteries/green-energy:** move to a single honest "**Current & forward interests**" line (not a claimed capability), or drop. Do not list capabilities without supporting work.
- Keep the writing factual and concrete; avoid filler superlatives.

### Task I — Add high-value content that is currently missing or buried
Confirm specifics with the owner / his CV, then add:
- **Open-source tools, front and center** (these matter a lot for an ML-for-science portfolio): `mlip_benchmark`, `MolSpecPy`, and contributions to the **hippynn** ecosystem. Give each a one-line description + GitHub link.
- **Conference leadership:** **co-organizer of MLCM-26** (Machine Learning in Chemical & Materials Sciences, Santa Fe) — a strong, real leadership signal.
- **Funded research:** the **$7.5M DOD-funded** collaborative project (polariton chemistry) he contributed to.
- **Education completeness:** add the **second M.Sc. (NIT Rourkela)** if not present.
- These should appear in the Projects/Open-Source, About, and Experience sections respectively.

### Task J — Tighten the "Technical Expertise" tool list to an honest, focused stack
- The current list reads like a kitchen sink (e.g. TensorFlow, ChemAxon, KNIME, SIESTA, NAMD). Trim to tools he genuinely uses day-to-day so the list is credible. Core, verifiable stack: **Python, C++, CUDA; PyTorch; HIP-NN/hippynn; RDKit, OpenBabel; TeraChem, Gaussian, ORCA, AMBER, OpenMM; SLURM, Ray; MOE (drug discovery).** Confirm the final list with the owner; remove anything he would not want to be interviewed on.

### Task K — Sync certifications
- Current certs (Udemy "Cheminformatics", Simplilearn "Data Science") read weak and may be dated. Confirm the current set with the owner and update (e.g. include the NVIDIA CUDA Python cert and any newer ML/AI engineering certificate he holds). Remove stale ones.

---

## 4. MINOR polish / housekeeping
- **SEO/meta:** ensure a good `<title>`, `<meta name="description">`, and **Open Graph / Twitter card** tags (title, description, image) so shared links preview well. Add a `favicon`.
- **Accessibility:** add descriptive `alt` text to all images (profile photo, project thumbnails); ensure heading levels are sequential; check color contrast on the colored stat tiles.
- **Responsive checks:** the timeline says "scroll horizontally to explore" — verify it works on mobile; test the whole page at narrow widths (the hero animation + stat grid especially).
- **Image hygiene:** confirm project images load (no 404s) and are reasonably sized/compressed.
- **Links audit:** click every external link (Scholar, ORCID, ResearchGate, GitHub, LinkedIn, X `@samdig`, project repos, tutorial repos) and fix any dead ones.
- **Performance:** defer non-critical JS; lazy-load below-the-fold images.
- **Consistency:** one consistent name for the CV/resume concept across nav, section, and CTA ("CV").

---

## 5. SUGGESTED ORDER & COMMITS
Work in small, reviewable commits on a feature branch (not directly on the default branch); open a PR or let the owner review before merging:
1. `fix: replace broken contact form with Google Form link` (Task A)
2. `fix: gate CV behind contact form, remove public resume` (Task B)
3. `fix: reconcile citation counts / live Scholar link` (Task C, D)
4. `content: reposition research focus to ML-for-science pillars` (Task F, H)
5. `content: add open-source tools, MLCM-26, DOD project, education` (Task I)
6. `content: trim tech stack + refresh certs` (Task J, K)
7. `chore: SEO/meta, a11y alt text, link audit, responsive fixes` (Task 4)
8. `chore: (optional) Tailwind build, page split` (Task E, G)
- After merge, verify the live site at https://ajaykhanna.github.io/ (GitHub Pages redeploys on push to the default branch).

## 6. ACCEPTANCE CHECKLIST (definition of done)
- [ ] Contact uses the owner's Google Form; no dead Formspree form remains.
- [ ] No public CV/resume is downloadable or directly reachable; CV is "available on request" via the form; `resume.html` deleted (or noindexed + unlinked, per owner).
- [ ] Citation numbers are consistent everywhere (or replaced by a live Scholar link); no "2000 downloads" filler; in-press paper labeled.
- [ ] Hero + research section lead with ML-for-science / MLIPs; battery/green-energy are not presented as established expertise.
- [ ] Open-source tools, MLCM-26, DOD project, full education present and accurate.
- [ ] Tech stack and certs trimmed to an honest, confirmed list.
- [ ] SEO meta + OG tags + favicon present; all images have alt text; all links work; layout holds on mobile.
- [ ] Owner has reviewed copy for factual accuracy before merge.

---

## 7. GROUND-TRUTH FACTS (use for accurate copy; confirm anything marked "verify")
> These are drawn from his current site and public profile. The owner is the source of truth — have him confirm details, especially counts and tool lists.

- **Name / role:** Ajay Khanna — Computational Chemist; Postdoctoral Researcher, Los Alamos National Laboratory (advisor: Dr. Sergei Tretiak), since Nov 2024.
- **Prior:** PhD, Computational Chemistry, UC Merced (Isborn group), 2024; M.Sc. UC Merced; **M.Sc., NIT Rourkela (verify)**; B.Sc. (Hons.) Chemistry, University of Delhi. CADD intern, Frontier Medicines (2023, BTK inhibitors).
- **Publications (from the site — verify counts via Scholar):**
  - "Covalent Control of Excitonic Interactions in Perylene Diimide Trimers," *Nano Letters*, 2026 (in press / accepted).
  - "Deconstructing Chirality: … Azobenzene Derivatives with X-ray Circular Dichroism," *J. Phys. Chem. Lett.*, 2025 (journal cover).
  - "Calculating Absorption and Fluorescence Spectra for Chromophores in Solution," *J. Chem. Phys.*, 2024.
  - "Axial H-Bonding Solvent Controls Inhomogeneous Spectral Broadening," *J. Phys. Chem. B*, 2024.
  - "Molecular Polariton Electroabsorption," *Nature Communications*, 13, 7937, 2022 (co-author).
  - "Explicit Environmental and Vibronic Effects in Simulations of … Optical Spectroscopy," *J. Chem. Phys.*, 2021.
  - Earlier: "Ligand Driven Electron Counting Rule Selection: Ge5R Complex," 2018.
- **Open-source tools:** `mlip_benchmark` (ML interatomic potential benchmarking), `MolSpecPy` (spectroscopy, TeraChem-interfaced), contributions to **hippynn**. (Plus existing site projects: QM/MM automation, MolVizMan, Excitonic Coupling, DOI2BibTex.)
- **Leadership / funding:** co-organizer, **MLCM-26** (Santa Fe); contributor to a **$7.5M DOD-funded** polariton-chemistry project. (verify amounts/labels with owner)
- **Core skills (honest stack):** ML interatomic potentials / HIP-NN(hippynn), GNN / message-passing, active learning, uncertainty quantification; QM/MM, DFT/TDDFT, (non)adiabatic molecular dynamics, spectroscopy; CADD (docking, MD, binding free energy, MOE); Python/C++/CUDA, PyTorch, RDKit, TeraChem, Gaussian, ORCA, AMBER, OpenMM, SLURM, Ray, GPU/HPC.
- **Profiles:** Google Scholar `user=qJM0sOIAAAAJ`, ORCID `0000-0002-8313-1393`, ResearchGate `Ajay-Khanna-2`, GitHub `Ajaykhanna`, LinkedIn `ajay-khanna`, X `@samdig`.
- **Contact form URL:** `{{GOOGLE_FORM_URL}}` — ask the owner.

---

*End of brief. Start with Section 1 (Tasks A–D). Ask the owner for the Google Form URL and the final tool/cert lists before writing the affected copy.*
