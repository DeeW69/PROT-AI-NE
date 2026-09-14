/**
 * PROT-AI-NE viewer (Couche 6).
 *
 * Consumes the JSON produced by `export.py` (protein_id + a list of
 * candidates, each with a rank, sequence, dot-bracket structure, MFE and
 * detailed scores). Two independent controls: a top1..topN selector and a
 * 2D/3D mode toggle.
 *
 * 2D mode renders the secondary structure (already computed, local, fast)
 * via fornac. 3D mode looks for a PDB file produced by the not-yet-built
 * `prediction_3d.py` (RhoFold+) at a conventional sibling path; per the
 * project's design this must degrade gracefully (never crash) when that
 * file is absent, unreachable, or the JSON was loaded from local disk
 * (file picker) rather than fetched from a URL.
 */
(function () {
  "use strict";

  const state = {
    data: null,
    baseUrl: null,
    rank: 1,
    mode: "2d",
    manualPdbText: null,
    manualPdbRank: null,
  };

  const els = {};

  function $(id) {
    return document.getElementById(id);
  }

  function init() {
    els.jsonFile = $("json-file-input");
    els.jsonUrl = $("json-url-input");
    els.jsonUrlBtn = $("json-url-load-btn");
    els.loadError = $("load-error");
    els.controls = $("controls");
    els.rankSelect = $("rank-select");
    els.proteinLabel = $("protein-label");
    els.mode2dBtn = $("mode-2d-btn");
    els.mode3dBtn = $("mode-3d-btn");
    els.render2d = $("render-2d");
    els.render3d = $("render-3d");
    els.render3dMessage = $("render-3d-message");
    els.pdbFileInput = $("pdb-file-input");
    els.scoreFitness = $("score-fitness");
    els.scoreGc = $("score-gc");
    els.scoreMfe = $("score-mfe");
    els.scoreCai = $("score-cai");
    els.fornaHost = els.render2d.querySelector(".forna-host");
    els.canvas3d = els.render3d.querySelector(".viewer-3d-canvas");

    els.jsonFile.addEventListener("change", onJsonFileChosen);
    els.jsonUrlBtn.addEventListener("click", onJsonUrlLoad);
    els.rankSelect.addEventListener("change", onRankChange);
    els.mode2dBtn.addEventListener("click", () => setMode("2d"));
    els.mode3dBtn.addEventListener("click", () => setMode("3d"));
    els.pdbFileInput.addEventListener("change", onPdbFileChosen);
  }

  // ---- loading -------------------------------------------------------

  function onJsonFileChosen(evt) {
    const file = evt.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      try {
        loadData(JSON.parse(reader.result), null);
      } catch (err) {
        showError("Fichier JSON invalide : " + err.message);
      }
    };
    reader.onerror = () => showError("Impossible de lire le fichier local.");
    reader.readAsText(file);
  }

  function onJsonUrlLoad() {
    const url = els.jsonUrl.value.trim();
    if (!url) return;
    fetch(url)
      .then((resp) => {
        if (!resp.ok) throw new Error("HTTP " + resp.status);
        return resp.json();
      })
      .then((data) => loadData(data, url))
      .catch((err) => {
        showError(
          "Impossible de charger '" + url + "' (" + err.message + "). " +
          "Rappel : ce chemin est relatif a viewer/index.html, pas au " +
          "dossier ou tourne le serveur - si le JSON est a la racine du " +
          "projet, il faut '../' devant (ex. '../top_candidates_structures.json'). " +
          "Si la page est ouverte directement depuis le disque (file://), " +
          "utilisez plutot le selecteur de fichier ci-dessus."
        );
      });
  }

  function showError(message) {
    els.loadError.textContent = message;
    els.loadError.hidden = !message;
  }

  function loadData(data, baseUrl) {
    if (!data || !Array.isArray(data.candidates) || data.candidates.length === 0) {
      showError("JSON charge mais aucun candidat trouve (structure inattendue).");
      return;
    }
    showError("");
    state.data = data;
    state.baseUrl = baseUrl;
    state.rank = data.candidates[0].rank;
    state.manualPdbText = null;
    state.manualPdbRank = null;

    els.proteinLabel.textContent = data.protein_id || "(proteine sans identifiant)";

    els.rankSelect.innerHTML = "";
    data.candidates.forEach((candidate) => {
      const option = document.createElement("option");
      option.value = String(candidate.rank);
      option.textContent = "Top " + candidate.rank;
      els.rankSelect.appendChild(option);
    });
    els.rankSelect.value = String(state.rank);

    els.controls.hidden = false;
    renderCurrent();
  }

  // ---- controls --------------------------------------------------------

  function onRankChange() {
    state.rank = parseInt(els.rankSelect.value, 10);
    renderCurrent();
  }

  function setMode(mode) {
    state.mode = mode;
    els.mode2dBtn.classList.toggle("active", mode === "2d");
    els.mode3dBtn.classList.toggle("active", mode === "3d");
    els.render2d.hidden = mode !== "2d";
    els.render3d.hidden = mode !== "3d";
    renderCurrent();
  }

  function currentCandidate() {
    if (!state.data) return null;
    return state.data.candidates.find((c) => c.rank === state.rank) || null;
  }

  // ---- scores ------------------------------------------------------

  function gcContent(sequence) {
    if (!sequence) return null;
    const upper = sequence.toUpperCase();
    let gc = 0;
    for (let i = 0; i < upper.length; i += 1) {
      if (upper[i] === "G" || upper[i] === "C") gc += 1;
    }
    return upper.length ? (100 * gc) / upper.length : null;
  }

  function fmt(value, digits, suffix) {
    if (value === null || value === undefined || Number.isNaN(value)) return "N/A";
    return value.toFixed(digits) + (suffix || "");
  }

  function renderScores(candidate) {
    const scores = candidate.scores || {};
    els.scoreFitness.textContent = fmt(scores.fitness, 4);
    els.scoreGc.textContent = fmt(gcContent(candidate.sequence), 1, " %");
    els.scoreMfe.textContent = fmt(candidate.mfe, 2, " kcal/mol");
    els.scoreCai.textContent = fmt(scores.cai, 4);
  }

  // ---- rendering ---------------------------------------------------

  function renderCurrent() {
    const candidate = currentCandidate();
    if (!candidate) return;
    renderScores(candidate);
    if (state.mode === "2d") {
      render2d(candidate);
    } else {
      render3d(candidate);
    }
  }

  function render2d(candidate) {
    els.fornaHost.innerHTML = "";
    const container = new fornac.FornaContainer(els.fornaHost, {
      applyForce: true,
      allowPanningAndZooming: true,
      labelInterval: 10,
    });
    container.addRNA(candidate.dot_bracket, {
      structure: candidate.dot_bracket,
      sequence: candidate.sequence,
    });
  }

  function pdbUrlFor(candidate) {
    if (!state.baseUrl) return null;
    // "<dir>/<stem>_structures.json" -> "<dir>/<stem>_structures3d/candidate_<rank>.pdb"
    // matching export.py's "<stem>_structures.json" naming convention.
    const match = state.baseUrl.match(/^(.*\/)?([^/]*?)(_structures)?\.json$/i);
    if (!match) return null;
    const dir = match[1] || "";
    const stem = match[2];
    return dir + stem + "_structures3d/candidate_" + candidate.rank + ".pdb";
  }

  function show3dUnavailable(message) {
    els.canvas3d.hidden = true;
    els.canvas3d.innerHTML = "";
    els.render3dMessage.hidden = false;
    els.render3dMessage.textContent = message;
  }

  function render3dFromText(pdbText) {
    els.render3dMessage.hidden = true;
    els.canvas3d.hidden = false;
    els.canvas3d.innerHTML = "";
    const viewer = $3Dmol.createViewer(els.canvas3d, { backgroundColor: "white" });
    viewer.addModel(pdbText, "pdb");
    viewer.setStyle({}, { cartoon: { color: "spectrum" } });
    viewer.zoomTo();
    viewer.render();
  }

  function render3d(candidate) {
    if (state.manualPdbText && state.manualPdbRank === candidate.rank) {
      render3dFromText(state.manualPdbText);
      return;
    }

    const url = pdbUrlFor(candidate);
    if (!url) {
      show3dUnavailable(
        "Structure 3D non disponible pour ce candidat (aucune prediction RhoFold+ " +
        "trouvee - la couche de prediction 3D n'est pas encore integree a ce run, " +
        "ou le JSON a ete charge depuis un fichier local). Vous pouvez charger " +
        "manuellement un fichier PDB ci-dessous si vous en avez genere un."
      );
      return;
    }

    fetch(url)
      .then((resp) => {
        if (!resp.ok) throw new Error("HTTP " + resp.status);
        return resp.text();
      })
      .then((text) => render3dFromText(text))
      .catch(() => {
        show3dUnavailable(
          "Structure 3D non disponible pour ce candidat (aucun fichier trouve a '" +
          url + "'). Vous pouvez charger manuellement un fichier PDB ci-dessous " +
          "si vous en avez genere un."
        );
      });
  }

  function onPdbFileChosen(evt) {
    const file = evt.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      state.manualPdbText = reader.result;
      state.manualPdbRank = state.rank;
      if (state.mode === "3d") {
        render3dFromText(state.manualPdbText);
      }
    };
    reader.readAsText(file);
  }

  document.addEventListener("DOMContentLoaded", init);
})();
