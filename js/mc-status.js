(async function () {
  const statusEl = document.getElementById("mc-status");
  if (!statusEl) return;

  const textEl = statusEl.querySelector(".status-text");

  try {
    const response = await fetch(
      "https://api.mcsrvstat.us/2/mc.cjones.dev"
    );
    const data = await response.json();

    if (data.online) {
      statusEl.classList.add("online");
      textEl.textContent = `Online — ${data.players.online}/${data.players.max} players · ${data.version}`;
    } else {
      statusEl.classList.add("offline");
      textEl.textContent = "Offline";
    }
  } catch (error) {
    statusEl.classList.add("offline");
    textEl.textContent = "Status unavailable";
  }
})();
