# Training — mes formations personnelles

Ce dépôt regroupe mes cours et supports de formation personnelle. Chaque cours vit
dans son propre dossier et peut se lire directement dans le navigateur.

## Cours disponibles

| Cours | Description | En ligne | Fichiers |
|-------|-------------|----------|----------|
| 🦀 **Rust pour la robotique temps réel** | Débuter Rust jusqu'aux implémentations de filtre de Kalman et de cinématique directe/inverse en temps réel. 14 chapitres en français : théorie, code commenté et exercices à solutions repliables. | [🌐 Voir le cours en ligne](https://core-engineering.github.io/training/cours-rust-robotique/index.html) | [`index.html`](cours-rust-robotique/index.html) |
| 🤖 **Cinématique robotique sur automate (PLC)** | Réviser la cinématique des bras manipulateurs — paramétrisation DH/MDH et autres, cinématique directe, Jacobien et calibration — avec une implémentation en IEC 61131-3 (Structured Text). 7 chapitres en français : théorie, code commenté et exercices à solutions repliables. | [🌐 Voir le cours en ligne](https://core-engineering.github.io/training/cours-robotique-cinematique/index.html) | [`index.html`](cours-robotique-cinematique/index.html) |
| 📈 **Ruckig — génération de trajectoire en ligne** | Comprendre l'algorithme de génération de trajectoire temps réel limitée en jerk : profil en S à 7 segments, familles de profil, solveur temps-optimal (Step 1 / Step 2), synchronisation multi-axes et re-ciblage à la volée. 8 chapitres en français, orientés théorie et illustrés par des graphiques cinématiquement exacts (sans exercices). Accompagne le portage SCL [`ruckig-scl`](https://github.com/pantor/ruckig). | [🌐 Voir le cours en ligne](https://core-engineering.github.io/training/cours-ruckig-otg/index.html) | [`index.html`](cours-ruckig-otg/index.html) |

*D'autres cours viendront s'ajouter ici au fil de mes formations.*

## Comment lire un cours

Ouvre le fichier `index.html` du dossier concerné dans ton navigateur (double-clic).
La coloration syntaxique du code s'appuie sur un CDN : avec une connexion internet
elle est plus jolie, mais le code reste parfaitement lisible hors-ligne.
