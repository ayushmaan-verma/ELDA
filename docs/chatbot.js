/**
 * chatbot.js
 * ─────────────────────────────────────────────────────────────────────────────
 * ELDA Documentation Chatbot — Core Engine
 *
 * Responsibilities:
 *  1. Build and tear down the chat widget DOM (with minimize/maximize/delete options)
 *  2. Save and restore conversation history via localStorage
 *  3. Accept user input (text + Enter key) and related question chips
 *  4. Run refined query matching with keyword weights & character fuzzy similarity
 *  5. Render bot replies with timestamps, code blocks (with Copy), warnings, and notes
 *  6. Manage custom modal confirmation dialog for deleting history
 *  7. Accessibility: Escape key, tab-navigation, and modal focus trapping
 *
 * Dependencies:  chatbot-data.js  (must be loaded first)
 * ─────────────────────────────────────────────────────────────────────────────
 */

(function () {
  "use strict";

  // ── Constants ──────────────────────────────────────────────────────────────

  /** Typing indicator delay in milliseconds (makes the bot feel more natural). */
  const TYPING_DELAY_MS = 600;

  /** Minimum fuzzy-match score (0–1) required to accept a result. */
  const FUZZY_THRESHOLD = 0.35;

  // ── State ──────────────────────────────────────────────────────────────────

  /** Whether the chatbot window is currently open. */
  let isOpen = false;

  /** Whether the welcome message has already been shown this session. */
  let welcomeShown = false;

  /** Whether the chatbot is currently minimized. */
  let isMinimized = false;

  /** Whether the chatbot is currently maximized. */
  let isMaximized = false;

  /** Stored array of messages for localStorage persistence. */
  let messageHistory = [];

  // ── DOM References (populated in buildWidget) ──────────────────────────────
  let widget          = null; // outer container #elda-chatbot
  let toggleBtn       = null; // floating toggle button
  let chatWindow      = null; // the sliding chat panel
  let messagesEl      = null; // scrollable messages area
  let inputEl         = null; // text input
  let sendBtn         = null; // send button
  let deleteBtn       = null; // delete history action button
  let minimizeBtn     = null; // minimize action button
  let maximizeBtn     = null; // maximize action button
  let confirmOverlay  = null; // custom confirmation dialog overlay

  // ── Initialisation ─────────────────────────────────────────────────────────

  /**
   * Entry point. Builds the widget and attaches it to <body>.
   * Called automatically when the script loads (see bottom of file).
   */
  function init() {
    buildWidget();
    attachEvents();
    loadHistory();
  }

  // ── Widget Construction ────────────────────────────────────────────────────

  /**
   * Creates all chatbot DOM nodes and appends them to document.body.
   * The widget starts hidden; toggling the button reveals it.
   */
  function buildWidget() {
    // ── Outer wrapper ──────────────────────────────────────────────────────
    widget = el("div", { id: "elda-chatbot", "aria-live": "polite" });

    // ── Floating toggle button ─────────────────────────────────────────────
    toggleBtn = el("button", {
      id: "elda-chat-toggle",
      type: "button",
      "aria-label": "Open ELDA documentation chatbot",
      "aria-expanded": "false",
      "aria-controls": "elda-chat-window",
      title: "Ask ELDA"
    });
    toggleBtn.innerHTML = `
      <span class="elda-chat-toggle-icon" aria-hidden="true">
        <!-- Chat bubble icon (open state) -->
        <svg class="icon-chat" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
             stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
          <path d="M21 15a2 2 0 0 1-2 2H7l-4 4V5a2 2 0 0 1 2-2h14a2 2 0 0 1 2 2z"/>
        </svg>
        <!-- X icon (close state) -->
        <svg class="icon-close" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
             stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round">
          <line x1="18" y1="6" x2="6" y2="18"/>
          <line x1="6" y1="6" x2="18" y2="18"/>
        </svg>
      </span>
      <span class="elda-chat-toggle-label">Ask ELDA</span>
    `;

    // ── Notification badge (pulse dot shown when closed) ──────────────────
    const badge = el("span", { class: "elda-chat-badge", "aria-hidden": "true" });
    toggleBtn.appendChild(badge);

    // ── Chat window ────────────────────────────────────────────────────────
    chatWindow = el("div", {
      id: "elda-chat-window",
      role: "dialog",
      "aria-modal": "false",
      "aria-label": "ELDA Documentation Chatbot"
    });

    // Header
    const header = el("header", { class: "elda-chat-header" });
    header.innerHTML = `
      <div class="elda-chat-header-info">
        <div class="elda-chat-avatar" aria-hidden="true">
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
               stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <rect x="3" y="3" width="18" height="18" rx="2"/>
            <path d="M9 9h6M9 12h6M9 15h4"/>
          </svg>
        </div>
        <div>
          <p class="elda-chat-name">ELDA Assistant</p>
          <p class="elda-chat-status">
            <span class="elda-status-dot" aria-hidden="true"></span>
            Always online
          </p>
        </div>
      </div>
      <div class="elda-chat-header-actions">
        <button class="elda-chat-action-btn elda-chat-delete" type="button" aria-label="Clear chat history" title="Clear Chat History">
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
               stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <polyline points="3 6 5 6 21 6"></polyline>
            <path d="M19 6v14a2 2 0 0 1-2 2H7a2 2 0 0 1-2-2V6m3 0V4a2 2 0 0 1 2-2h4a2 2 0 0 1 2 2v2"></path>
          </svg>
        </button>
        <button class="elda-chat-action-btn elda-chat-minimize" type="button" aria-label="Minimize chatbot" title="Minimize">
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
               stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <line x1="5" y1="12" x2="19" y2="12"></line>
          </svg>
        </button>
        <button class="elda-chat-action-btn elda-chat-maximize" type="button" aria-label="Maximize chatbot" title="Maximize">
          <svg class="icon-max" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
               stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <rect x="3" y="3" width="18" height="18" rx="2" ry="2"></rect>
          </svg>
          <svg class="icon-restore" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
               stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <path d="M4 14h6v6m10-10h-6V4"></path>
          </svg>
        </button>
        <button class="elda-chat-action-btn elda-chat-close" type="button" aria-label="Close chatbot" title="Close">
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
               stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round">
            <line x1="18" y1="6" x2="6" y2="18"/>
            <line x1="6" y1="6" x2="18" y2="18"/>
          </svg>
        </button>
      </div>
    `;

    deleteBtn   = header.querySelector(".elda-chat-delete");
    minimizeBtn = header.querySelector(".elda-chat-minimize");
    maximizeBtn = header.querySelector(".elda-chat-maximize");

    // Messages area
    messagesEl = el("div", {
      class: "elda-chat-messages",
      id: "elda-chat-messages",
      role: "log",
      "aria-label": "Chat messages"
    });

    // Quick-start chips (shown once at top)
    const quickStart = buildQuickStart();

    // Input area
    const inputArea = el("div", { class: "elda-chat-input-area" });
    inputEl = el("input", {
      type: "text",
      class: "elda-chat-input",
      id: "elda-chat-input",
      placeholder: "Ask about ELDA…",
      autocomplete: "off",
      "aria-label": "Type your question",
      maxlength: "300"
    });
    sendBtn = el("button", {
      type: "button",
      class: "elda-chat-send",
      "aria-label": "Send message",
      title: "Send"
    });
    sendBtn.innerHTML = `
      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
           stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
        <line x1="22" y1="2" x2="11" y2="13"/>
        <polygon points="22 2 15 22 11 13 2 9 22 2"/>
      </svg>
    `;

    inputArea.append(inputEl, sendBtn);

    // Footer branding
    const footer = el("footer", { class: "elda-chat-footer" });
    footer.innerHTML = `<p>Powered by <strong>ELDA Docs</strong> · Offline-capable</p>`;

    // Custom Confirmation Dialog Overlay
    confirmOverlay = el("div", {
      class: "elda-confirm-overlay",
      id: "elda-confirm-overlay",
      role: "alertdialog",
      "aria-modal": "true",
      "aria-hidden": "true",
      "aria-describedby": "elda-confirm-msg"
    });
    confirmOverlay.innerHTML = `
      <div class="elda-confirm-modal">
        <p id="elda-confirm-msg">Are you sure you want to clear this conversation?</p>
        <div class="elda-confirm-buttons">
          <button class="elda-confirm-btn elda-confirm-yes" type="button" tabindex="0">Yes, Clear</button>
          <button class="elda-confirm-btn elda-confirm-cancel" type="button" tabindex="0">Cancel</button>
        </div>
      </div>
    `;

    // Assemble chat window
    chatWindow.append(header, messagesEl, quickStart, inputArea, footer, confirmOverlay);

    // Assemble widget
    widget.append(toggleBtn, chatWindow);
    document.body.appendChild(widget);
  }

  /**
   * Builds the quick-start chip strip shown below the welcome message.
   * @returns {HTMLElement}
   */
  function buildQuickStart() {
    const strip = el("div", { class: "elda-chat-quickstart", id: "elda-chat-quickstart" });
    const label = el("p", { class: "elda-quickstart-label" });
    label.textContent = "Suggested questions:";
    strip.appendChild(label);

    const chips = el("div", { class: "elda-quickstart-chips" });
    (typeof STARTER_QUESTIONS !== "undefined" ? STARTER_QUESTIONS : []).forEach((q) => {
      const chip = el("button", { type: "button", class: "elda-quickstart-chip", tabindex: "0" });
      chip.textContent = q;
      chip.addEventListener("click", () => {
        hideQuickStart();
        submitUserMessage(q);
      });
      chips.appendChild(chip);
    });

    strip.appendChild(chips);
    return strip;
  }

  // ── Event Wiring ──────────────────────────────────────────────────────────

  /** Attaches all interactive event listeners after the DOM is built. */
  function attachEvents() {
    // Toggle button opens/closes the chat window
    toggleBtn.addEventListener("click", toggleChat);

    // Close button inside the chat header actions
    chatWindow.querySelector(".elda-chat-close").addEventListener("click", (e) => {
      e.stopPropagation();
      closeChat();
    });

    // Minimize & Maximize action buttons
    minimizeBtn.addEventListener("click", toggleMinimize);
    maximizeBtn.addEventListener("click", toggleMaximize);

    // Delete chat history action button triggers custom confirmation
    deleteBtn.addEventListener("click", showConfirmDelete);

    // Confirmation dialog buttons
    confirmOverlay.querySelector(".elda-confirm-yes").addEventListener("click", handleConfirmClear);
    confirmOverlay.querySelector(".elda-confirm-cancel").addEventListener("click", hideConfirmDelete);

    // Header click restores the window if minimized
    const header = chatWindow.querySelector(".elda-chat-header");
    header.addEventListener("click", (e) => {
      if (isMinimized && !e.target.closest(".elda-chat-action-btn")) {
        restoreChat();
      }
    });

    // Send button click
    sendBtn.addEventListener("click", handleSend);

    // Enter key in input submits the message
    inputEl.addEventListener("keydown", (e) => {
      if (e.key === "Enter" && !e.shiftKey) {
        e.preventDefault();
        handleSend();
      }
    });

    // Show send button only when input has text
    inputEl.addEventListener("input", () => {
      sendBtn.classList.toggle("elda-send-active", inputEl.value.trim().length > 0);
    });

    // Close when clicking the dark overlay/background (mobile focus loss)
    document.addEventListener("click", (e) => {
      if (
        isOpen &&
        !chatWindow.contains(e.target) &&
        !toggleBtn.contains(e.target)
      ) {
        closeChat();
      }
    });

    // Accessible keyboard navigation - Global Keydowns
    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape") {
        if (confirmOverlay.classList.contains("elda-confirm-show")) {
          hideConfirmDelete();
        } else if (isOpen) {
          closeChat();
        }
      }
    });

    // Modal focus trapping inside the confirmation overlay
    confirmOverlay.addEventListener("keydown", (e) => {
      if (e.key === "Tab") {
        const focusables = confirmOverlay.querySelectorAll("button");
        const first = focusables[0];
        const last = focusables[focusables.length - 1];
        if (e.shiftKey) {
          if (document.activeElement === first) {
            last.focus();
            e.preventDefault();
          }
        } else {
          if (document.activeElement === last) {
            first.focus();
            e.preventDefault();
          }
        }
      }
    });
  }

  // ── Open / Close / Minimize / Maximize ────────────────────────────────────

  /** Toggles the chat window open or closed. */
  function toggleChat() {
    isOpen ? closeChat() : openChat();
  }

  /** Opens the chat window and shows the welcome message on first open. */
  function openChat() {
    isOpen = true;
    chatWindow.classList.add("elda-chat-open");
    toggleBtn.setAttribute("aria-expanded", "true");
    toggleBtn.classList.add("elda-toggle-open");
    // Remove the notification badge once opened
    toggleBtn.querySelector(".elda-chat-badge").classList.remove("elda-badge-pulse");

    // Show welcome message only once (if history is empty)
    if (!welcomeShown && messageHistory.length === 0) {
      welcomeShown = true;
      showTypingIndicator();
      setTimeout(() => {
        removeTypingIndicator();
        addBotMessage(
          "Hello! 👋 I'm ELDA Assistant. I can help you explore the ELDA documentation, explain matrix operations, guide you through installation, and answer questions about the library. Try asking: 'How do I build ELDA?' or 'What is matrix.hpp?'"
        );
      }, 300);
    }

    // Focus input for accessibility
    setTimeout(() => {
      if (!isMinimized) inputEl.focus();
    }, 300);
  }

  /** Closes the chat window. */
  function closeChat() {
    isOpen = false;
    chatWindow.classList.remove("elda-chat-open");
    toggleBtn.setAttribute("aria-expanded", "false");
    toggleBtn.classList.remove("elda-toggle-open");
    // If confirmation is open when closed, hide it
    hideConfirmDelete();
  }

  /** Toggles the minimized state of the chat window. */
  function toggleMinimize(e) {
    if (e) e.stopPropagation();
    isMinimized ? restoreChat() : minimizeChat();
  }

  /** Collapses the chat window to show only the header bar. */
  function minimizeChat() {
    isMinimized = true;
    chatWindow.classList.add("elda-chat-minimized");
    chatWindow.classList.remove("elda-chat-maximized");
    isMaximized = false;
    minimizeBtn.setAttribute("aria-label", "Restore chatbot");
    minimizeBtn.setAttribute("title", "Restore");
    // Change minimize icon to a toggle Chevron Up
    minimizeBtn.innerHTML = `
      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
           stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
        <polyline points="18 15 12 9 6 15"></polyline>
      </svg>
    `;
    inputEl.blur();
  }

  /** Restores the minimized chat window back to normal size. */
  function restoreChat() {
    isMinimized = false;
    chatWindow.classList.remove("elda-chat-minimized");
    minimizeBtn.setAttribute("aria-label", "Minimize chatbot");
    minimizeBtn.setAttribute("title", "Minimize");
    // Restore minimize icon to a standard minus/line
    minimizeBtn.innerHTML = `
      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
           stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
        <line x1="5" y1="12" x2="19" y2="12"></line>
      </svg>
    `;
    inputEl.focus();
  }

  /** Toggles the maximized state of the chat window. */
  function toggleMaximize(e) {
    if (e) e.stopPropagation();
    if (isMaximized) {
      isMaximized = false;
      chatWindow.classList.remove("elda-chat-maximized");
      maximizeBtn.setAttribute("aria-label", "Maximize chatbot");
      maximizeBtn.setAttribute("title", "Maximize");
      chatWindow.classList.remove("elda-max-mode");
    } else {
      isMaximized = true;
      isMinimized = false;
      chatWindow.classList.remove("elda-chat-minimized");
      chatWindow.classList.add("elda-chat-maximized");
      chatWindow.classList.add("elda-max-mode");
      maximizeBtn.setAttribute("aria-label", "Restore chatbot size");
      maximizeBtn.setAttribute("title", "Restore Size");
      // Restore minimize button icon state
      minimizeBtn.innerHTML = `
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
             stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
          <line x1="5" y1="12" x2="19" y2="12"></line>
        </svg>
      `;
    }
  }

  /** Hides the quick-start chip strip (called after first interaction). */
  function hideQuickStart() {
    const qs = document.getElementById("elda-chat-quickstart");
    if (qs) {
      qs.style.maxHeight = "0";
      qs.style.opacity = "0";
      qs.style.pointerEvents = "none";
      setTimeout(() => qs.remove(), 300);
    }
  }

  // ── Message Handling & History ────────────────────────────────────────────

  /** Gets current time formatted as e.g. "1:45 PM". */
  function getFormattedTime() {
    const now = new Date();
    let hours = now.getHours();
    const minutes = String(now.getMinutes()).padStart(2, "0");
    const ampm = hours >= 12 ? "PM" : "AM";
    hours = hours % 12;
    hours = hours ? hours : 12; // 0 hour should be 12
    return `${hours}:${minutes} ${ampm}`;
  }

  /** Saves a message to the internal array and local storage. */
  function saveToHistory(sender, text, entryId = null, timestamp = null) {
    const time = timestamp || getFormattedTime();
    messageHistory.push({ sender, text, entryId, timestamp: time });
    try {
      localStorage.setItem("elda_chat_history", JSON.stringify(messageHistory));
    } catch (e) {
      console.warn("localStorage is not accessible:", e);
    }
  }

  /** Restores previous conversation from local storage. */
  function loadHistory() {
    try {
      const stored = localStorage.getItem("elda_chat_history");
      if (stored) {
        messageHistory = JSON.parse(stored);
        if (messageHistory.length > 0) {
          welcomeShown = true;
          hideQuickStart();

          messageHistory.forEach((msg) => {
            if (msg.sender === "user") {
              renderUserMessageDom(msg.text, msg.timestamp);
            } else {
              if (msg.entryId) {
                const entry = CHATBOT_DATA.find((e) => e.id === msg.entryId);
                if (entry) {
                  renderBotAnswerDom(entry, msg.timestamp);
                } else {
                  renderBotMessageDom(msg.text, msg.timestamp);
                }
              } else {
                renderBotMessageDom(msg.text, msg.timestamp);
              }
            }
          });
          scrollToBottom();
        }
      }
    } catch (e) {
      console.warn("Failed to load chat history:", e);
    }
  }

  /** Triggers the confirmation overlay to delete history. */
  function showConfirmDelete(e) {
    if (e) e.stopPropagation();
    confirmOverlay.setAttribute("aria-hidden", "false");
    confirmOverlay.classList.add("elda-confirm-show");
    // Shift focus to confirmation Cancel for safety
    setTimeout(() => {
      confirmOverlay.querySelector(".elda-confirm-cancel").focus();
    }, 100);
  }

  /** Dismisses the confirmation overlay. */
  function hideConfirmDelete() {
    confirmOverlay.setAttribute("aria-hidden", "true");
    confirmOverlay.classList.remove("elda-confirm-show");
    inputEl.focus();
  }

  /** Clear logic called after clicking Yes in confirmation dialog. */
  function handleConfirmClear() {
    confirmOverlay.setAttribute("aria-hidden", "true");
    confirmOverlay.classList.remove("elda-confirm-show");
    clearChatHistory();
  }

  /** Clears local storage history and displays success message. */
  function clearChatHistory() {
    messageHistory = [];
    try {
      localStorage.removeItem("elda_chat_history");
    } catch (e) {
      console.warn("Failed to clear local storage:", e);
    }

    // Clear messages UI
    messagesEl.innerHTML = "";

    // Requirement: On confirmation: "Chat history cleared successfully."
    addBotMessage("Chat history cleared successfully.");

    // Prompt welcome message again after a delay to keep the bot ready
    setTimeout(() => {
      welcomeShown = true;
      addBotMessage(
        "Hello! 👋 I'm ELDA Assistant. I can help you explore the ELDA documentation, explain matrix operations, guide you through installation, and answer questions about the library. Try asking: 'How do I build ELDA?' or 'What is matrix.hpp?'"
      );
    }, 1000);
  }

  /** Called when the user clicks Send or presses Enter. */
  function handleSend() {
    const text = inputEl.value.trim();
    if (!text) return;
    inputEl.value = "";
    sendBtn.classList.remove("elda-send-active");
    submitUserMessage(text);
  }

  /**
   * Appends the user's message bubble, then triggers the bot response.
   * @param {string} text - The user's message text.
   */
  function submitUserMessage(text) {
    hideQuickStart();
    addUserMessage(text);
    showTypingIndicator();

    // Small delay to simulate bot "thinking"
    setTimeout(() => {
      removeTypingIndicator();
      const entry = findAnswer(text);
      if (entry) {
        renderBotAnswer(entry);
      } else {
        renderFallbackReply();
      }
    }, TYPING_DELAY_MS);
  }

  // ── Matching Engine ───────────────────────────────────────────────────────

  /**
   * Main answer-lookup function. Compiles a composite score for each FAQ entry:
   *  1. Exact keyword phrase match (preserving stop words for greetings & short queries)
   *  2. Partial/contains match against keywords
   *  3. Token-based overlap ratio (useful keywords only)
   *  4. Character fuzzy similarity using Sørensen–Dice bigrams
   *
   * @param {string} userText - Raw user input.
   * @returns {Object|null} A CHATBOT_DATA entry or null if nothing qualifies.
   */
  function findAnswer(userText) {
    const rawQuery = userText
      .toLowerCase()
      .replace(/[^a-z0-9\s_]/g, " ")
      .replace(/\s+/g, " ")
      .trim();

    const query = normalize(userText);
    const queryTokens = query.split(/\s+/).filter(Boolean);

    let bestEntry = null;
    let highestScore = 0;

    for (const entry of CHATBOT_DATA) {
      let entryBestScore = 0;

      // ── Check all keywords for this entry ──
      for (const kw of entry.keywords) {
        const cleanKw = kw.toLowerCase().trim();
        let score = 0;

        // 1. Exact phrase match (weight: 10)
        if (rawQuery === cleanKw) {
          score = 10.0;
        }
        // 2. Contains match (weight: 5)
        else if (rawQuery.includes(cleanKw) || cleanKw.includes(rawQuery)) {
          score = 5.0;
        }
        // 3. Token-based overlap (weight: 3 * fraction matched)
        else if (queryTokens.length > 0) {
          const kwTokens = cleanKw.split(/\s+/).filter(Boolean);
          const matched = kwTokens.filter((t) =>
            queryTokens.some((qt) => qt.includes(t) || t.includes(qt))
          ).length;
          if (matched > 0) {
            score = 3.0 * (matched / kwTokens.length);
          }
        }

        // 4. Fuzzy bigram similarity (weight: 2 * similarity)
        const fuzzy = fuzzyScore_bigramSimilarity(rawQuery, cleanKw);
        if (fuzzy >= FUZZY_THRESHOLD) {
          const fuzzyScore = 2.0 * fuzzy;
          if (fuzzyScore > score) {
            score = fuzzyScore;
          }
        }

        if (score > entryBestScore) {
          entryBestScore = score;
        }
      }

      // ── Check canonical question text ──
      const cleanQuestion = entry.question
        .toLowerCase()
        .replace(/[^a-z0-9\s_]/g, " ")
        .replace(/\s+/g, " ")
        .trim();

      if (rawQuery === cleanQuestion) {
        entryBestScore = Math.max(entryBestScore, 10.0);
      } else if (rawQuery.includes(cleanQuestion) || cleanQuestion.includes(rawQuery)) {
        entryBestScore = Math.max(entryBestScore, 5.0);
      } else {
        const fuzzy = fuzzyScore_bigramSimilarity(rawQuery, cleanQuestion);
        if (fuzzy >= FUZZY_THRESHOLD) {
          entryBestScore = Math.max(entryBestScore, 2.0 * fuzzy);
        }
      }

      // Rank best overall matches
      if (entryBestScore > highestScore) {
        highestScore = entryBestScore;
        bestEntry = entry;
      }
    }

    const MIN_ACCEPTABLE_SCORE = 1.0;
    return highestScore >= MIN_ACCEPTABLE_SCORE ? bestEntry : null;
  }

  /**
   * Normalises a string for matching: lower-case, collapse whitespace,
   * strip common stop-words and punctuation.
   * @param {string} str
   * @returns {string}
   */
  function normalize(str) {
    const stopWords = new Set([
      "a","an","the","is","are","was","were","i","do","does","did","how",
      "what","where","when","why","which","who","can","could","would","should",
      "may","might","me","my","we","our","you","your","it","its","this",
      "that","these","those","in","on","at","to","for","of","and","or","but",
      "not","with","about","from","by","as","if","be","will","have","has","had"
    ]);

    return str
      .toLowerCase()
      .replace(/[^a-z0-9\s_]/g, " ") // strip punctuation
      .replace(/\s+/g, " ")
      .trim()
      .split(" ")
      .filter((w) => w.length > 1 && !stopWords.has(w))
      .join(" ");
  }

  /**
   * Computes bigram-based Sørensen–Dice similarity between two strings.
   * @param {string} a
   * @param {string} b
   * @returns {number}
   */
  function fuzzyScore_bigramSimilarity(a, b) {
    if (!a || !b) return 0;
    if (a === b) return 1;

    function bigrams(s) {
      const result = [];
      for (let i = 0; i < s.length - 1; i++) {
        result.push(s.slice(i, i + 2));
      }
      return result;
    }

    const aBigrams = bigrams(a);
    const bBigrams = bigrams(b);

    if (aBigrams.length === 0 || bBigrams.length === 0) return 0;

    const bSet = [...bBigrams];
    let intersection = 0;
    for (const bg of aBigrams) {
      const idx = bSet.indexOf(bg);
      if (idx !== -1) {
        intersection++;
        bSet.splice(idx, 1);
      }
    }

    return (2 * intersection) / (aBigrams.length + bBigrams.length);
  }

  // ── Rendering ─────────────────────────────────────────────────────────────

  /**
   * Appends a user message bubble to the chat and saves it to history.
   * @param {string} text
   */
  function addUserMessage(text) {
    const timestamp = getFormattedTime();
    renderUserMessageDom(text, timestamp);
    saveToHistory("user", text, null, timestamp);
  }

  /** Helper to append user message bubble to the DOM. */
  function renderUserMessageDom(text, timestamp) {
    const row = el("div", { class: "elda-msg-row elda-msg-user" });
    const bubble = el("div", { class: "elda-msg-bubble elda-msg-bubble-user" });

    const textSpan = el("span", { class: "elda-msg-text" });
    textSpan.textContent = text;

    const timeSpan = el("span", { class: "elda-msg-time" });
    timeSpan.textContent = timestamp;

    bubble.append(textSpan, timeSpan);
    row.appendChild(bubble);
    messagesEl.appendChild(row);
    scrollToBottom();
  }

  /**
   * Appends a plain bot text message bubble and saves it to history.
   * @param {string} text
   */
  function addBotMessage(text) {
    const timestamp = getFormattedTime();
    renderBotMessageDom(text, timestamp);
    saveToHistory("bot", text, null, timestamp);
  }

  /** Helper to append standard text bot bubble to the DOM. */
  function renderBotMessageDom(text, timestamp) {
    const row = el("div", { class: "elda-msg-row elda-msg-bot" });
    const avatar = el("div", { class: "elda-msg-avatar", "aria-hidden": "true" });
    avatar.innerHTML = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
      stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
      <rect x="3" y="3" width="18" height="18" rx="2"/>
      <path d="M9 9h6M9 12h6M9 15h4"/>
    </svg>`;

    const bubble = el("div", { class: "elda-msg-bubble elda-msg-bubble-bot" });

    const textDiv = el("div", { class: "elda-msg-text" });
    textDiv.innerHTML = renderMarkdownLite(text);

    const timeSpan = el("span", { class: "elda-msg-time" });
    timeSpan.textContent = timestamp;

    bubble.append(textDiv, timeSpan);
    row.append(avatar, bubble);
    messagesEl.appendChild(row);
    scrollToBottom();
  }

  /**
   * Renders a full FAQ entry as a bot message with optional components,
   * then pushes it to history.
   * @param {Object} entry - A FAQ database entry.
   */
  function renderBotAnswer(entry) {
    const timestamp = getFormattedTime();
    renderBotAnswerDom(entry, timestamp);
    saveToHistory("bot", entry.answer, entry.id, timestamp);
  }

  /** Helper to render structured bot reply blocks to the DOM. */
  function renderBotAnswerDom(entry, timestamp) {
    const row = el("div", { class: "elda-msg-row elda-msg-bot" });
    const avatar = el("div", { class: "elda-msg-avatar", "aria-hidden": "true" });
    avatar.innerHTML = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
      stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
      <rect x="3" y="3" width="18" height="18" rx="2"/>
      <path d="M9 9h6M9 12h6M9 15h4"/>
    </svg>`;

    const bubble = el("div", { class: "elda-msg-bubble elda-msg-bubble-bot" });

    // ── Main answer text ──
    const answerEl = el("div", { class: "elda-answer-text" });
    answerEl.innerHTML = renderMarkdownLite(entry.answer);
    bubble.appendChild(answerEl);

    // ── Warning badge ──
    if (entry.warning) {
      const warn = el("div", { class: "elda-badge elda-badge-warning", role: "alert" });
      warn.innerHTML = `<span class="elda-badge-icon" aria-hidden="true">⚠️</span>
                        <span>${escapeHtml(entry.warning)}</span>`;
      bubble.appendChild(warn);
    }

    // ── Note badge ──
    if (entry.note) {
      const note = el("div", { class: "elda-badge elda-badge-note" });
      note.innerHTML = `<span class="elda-badge-icon" aria-hidden="true">💡</span>
                        <span>${escapeHtml(entry.note)}</span>`;
      bubble.appendChild(note);
    }

    // ── Code block with Copy button ──
    if (entry.code) {
      const codeWrap = el("div", { class: "elda-code-block" });
      const copyBtn = el("button", {
        type: "button",
        class: "elda-code-copy",
        title: "Copy code",
        "aria-label": "Copy code snippet",
        tabindex: "0"
      });
      copyBtn.textContent = "Copy";
      copyBtn.addEventListener("click", () => copyCode(entry.code, copyBtn));

      const pre = el("pre", {});
      const code = el("code", {});
      code.textContent = entry.code;
      pre.appendChild(code);

      codeWrap.append(copyBtn, pre);
      bubble.appendChild(codeWrap);
    }

    // ── Related questions chips ──
    if (entry.related && entry.related.length > 0) {
      const relatedWrap = el("div", { class: "elda-related" });
      const relatedLabel = el("p", { class: "elda-related-label" });
      relatedLabel.textContent = "Related:";
      relatedWrap.appendChild(relatedLabel);

      const chips = el("div", { class: "elda-related-chips" });
      entry.related.forEach((relId) => {
        const relEntry = CHATBOT_DATA.find((e) => e.id === relId);
        if (!relEntry) return;
        const chip = el("button", { type: "button", class: "elda-related-chip", tabindex: "0" });
        chip.textContent = relEntry.question;
        chip.addEventListener("click", () => {
          addUserMessage(relEntry.question);
          showTypingIndicator();
          setTimeout(() => {
            removeTypingIndicator();
            renderBotAnswer(relEntry);
          }, TYPING_DELAY_MS);
        });
        chips.appendChild(chip);
      });

      relatedWrap.appendChild(chips);
      bubble.appendChild(relatedWrap);
    }

    // ── Timestamp ──
    const timeSpan = el("span", { class: "elda-msg-time" });
    timeSpan.textContent = timestamp;
    bubble.appendChild(timeSpan);

    row.append(avatar, bubble);
    messagesEl.appendChild(row);
    scrollToBottom();
  }

  /** Renders the fallback reply including 3 automatic suggestion chips. */
  function renderFallbackReply() {
    const timestamp = getFormattedTime();

    // Get 3 suggestions from starter questions
    const list = typeof STARTER_QUESTIONS !== "undefined" ? STARTER_QUESTIONS : [
      "What is ELDA?",
      "How do I build ELDA?",
      "How can I contribute?"
    ];
    const suggestions = [...list].sort(() => 0.5 - Math.random()).slice(0, 3);

    const row = el("div", { class: "elda-msg-row elda-msg-bot" });
    const avatar = el("div", { class: "elda-msg-avatar", "aria-hidden": "true" });
    avatar.innerHTML = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
      stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
      <rect x="3" y="3" width="18" height="18" rx="2"/>
      <path d="M9 9h6M9 12h6M9 15h4"/>
    </svg>`;

    const bubble = el("div", { class: "elda-msg-bubble elda-msg-bubble-bot" });

    // Main fallback message lines
    const textSpan = el("div", { class: "elda-msg-text" });
    textSpan.innerHTML = `Sorry, I couldn't find information related to that.<br>Please explore the documentation menu for more details.`;
    bubble.appendChild(textSpan);

    // Clickable suggested question chips
    const relatedWrap = el("div", { class: "elda-related" });
    const relatedLabel = el("p", { class: "elda-related-label" });
    relatedLabel.textContent = "Try asking:";
    relatedWrap.appendChild(relatedLabel);

    const chips = el("div", { class: "elda-related-chips" });
    suggestions.forEach((qText) => {
      const chip = el("button", { type: "button", class: "elda-related-chip", tabindex: "0" });
      chip.textContent = qText;
      chip.addEventListener("click", () => {
        submitUserMessage(qText);
      });
      chips.appendChild(chip);
    });
    relatedWrap.appendChild(chips);
    bubble.appendChild(relatedWrap);

    // Timestamp
    const timeSpan = el("span", { class: "elda-msg-time" });
    timeSpan.textContent = timestamp;
    bubble.appendChild(timeSpan);

    row.append(avatar, bubble);
    messagesEl.appendChild(row);
    scrollToBottom();

    // Save fallback state to history
    saveToHistory(
      "bot",
      "Sorry, I couldn't find information related to that. Please explore the documentation menu for more details.",
      null,
      timestamp
    );
  }

  // ── Typing Indicator ──────────────────────────────────────────────────────

  /** Shows the animated "bot is typing" indicator bubble. */
  function showTypingIndicator() {
    removeTypingIndicator(); // safety: ensure only one at a time
    const row = el("div", {
      class: "elda-msg-row elda-msg-bot elda-typing-row",
      id: "elda-typing-indicator"
    });
    const avatar = el("div", { class: "elda-msg-avatar", "aria-hidden": "true" });
    avatar.innerHTML = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
      stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
      <rect x="3" y="3" width="18" height="18" rx="2"/>
      <path d="M9 9h6M9 12h6M9 15h4"/>
    </svg>`;
    const bubble = el("div", { class: "elda-msg-bubble elda-msg-bubble-bot elda-typing-bubble" });
    bubble.setAttribute("aria-label", "ELDA bot is typing");
    bubble.innerHTML = `
      <span class="elda-typing-dots" aria-hidden="true">
        <span></span><span></span><span></span>
      </span>
    `;
    row.append(avatar, bubble);
    messagesEl.appendChild(row);
    scrollToBottom();
  }

  /** Removes the typing indicator from the DOM if present. */
  function removeTypingIndicator() {
    const indicator = document.getElementById("elda-typing-indicator");
    if (indicator) indicator.remove();
  }

  // ── Utilities ──────────────────────────────────────────────────────────────

  /**
   * Smoothly scrolls the messages container to the bottom.
   */
  function scrollToBottom() {
    requestAnimationFrame(() => {
      messagesEl.scrollTop = messagesEl.scrollHeight;
    });
  }

  /**
   * Copies a code string to the clipboard and gives visual feedback.
   */
  function copyCode(code, btn) {
    navigator.clipboard
      .writeText(code)
      .then(() => {
        btn.textContent = "Copied!";
        btn.classList.add("elda-code-copied");
        setTimeout(() => {
          btn.textContent = "Copy";
          btn.classList.remove("elda-code-copied");
        }, 1800);
      })
      .catch(() => {
        // Fallback for older browsers / non-HTTPS local context
        const ta = document.createElement("textarea");
        ta.value = code;
        ta.style.position = "fixed";
        ta.style.opacity = "0";
        document.body.appendChild(ta);
        ta.select();
        document.execCommand("copy");
        document.body.removeChild(ta);
        btn.textContent = "Copied!";
        setTimeout(() => (btn.textContent = "Copy"), 1800);
      });
  }

  /**
   * Renders a small subset of markdown:
   *  **bold**  →  <strong>bold</strong>
   *  `code`    →  <code>code</code>
   *  \n        →  <br>
   *  • lines   →  styled list items
   */
  function renderMarkdownLite(text) {
    if (!text) return "";

    let html = escapeHtml(text);

    // **bold**
    html = html.replace(/\*\*(.+?)\*\*/g, "<strong>$1</strong>");

    // `inline code`
    html = html.replace(/`([^`]+)`/g, "<code>$1</code>");

    // Bullet lines starting with • or -
    html = html.replace(/^[•\-]\s(.+)$/gm, '<span class="elda-list-item">$1</span>');

    // Newlines → <br>
    html = html.replace(/\n/g, "<br>");

    return html;
  }

  /**
   * Escapes HTML special characters to prevent XSS injection.
   */
  function escapeHtml(str) {
    return String(str)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

  /**
   * Convenience factory for creating DOM elements with attributes.
   */
  function el(tag, attrs = {}) {
    const node = document.createElement(tag);
    for (const [key, val] of Object.entries(attrs)) {
      node.setAttribute(key, val);
    }
    return node;
  }

  // ── Boot ──────────────────────────────────────────────────────────────────

  // Initialise after the DOM is ready
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }

  // Show the notification pulse on the toggle button after a short delay
  setTimeout(() => {
    if (!isOpen && toggleBtn) {
      const badgeEl = toggleBtn.querySelector(".elda-chat-badge");
      if (badgeEl && !localStorage.getItem("elda_chat_history")) {
        badgeEl.classList.add("elda-badge-pulse");
      }
    }
  }, 2000);
})();
