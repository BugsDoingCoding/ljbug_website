const SERVER_ADDRESS = "mc.cjones.dev";
const statusEl = document.getElementById("mc-status");

fetch(`https://api.mcstatus.io/v2/status/java/${SERVER_ADDRESS}`)
  .then(res => res.json())
  .then(data => {
    if (!statusEl) return;

    if (data.online) {
      const players = data.players.online;
      const max = data.players.max;

      statusEl.textContent = `🟢 Online — ${players}/${max} players`;
      statusEl.classList.remove("loading");
      statusEl.classList.add("online");
    } else {
      statusEl.textContent = "🔴 Offline";
      statusEl.classList.remove("loading");
      statusEl.classList.add("offline");
    }
  })
  .catch(() => {
    statusEl.textContent = "⚠️ Unable to fetch server status";
    statusEl.classList.remove("loading");
    statusEl.classList.add("offline");
  });
