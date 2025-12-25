const SERVER = "mc.cjones.dev";
const statusEl = document.getElementById("mc-status");
const motdEl = document.getElementById("mc-motd");

async function updateStatus() {
  try {
    const res = await fetch(`https://api.mcstatus.io/v2/status/java/${SERVER}`);
    const data = await res.json();

    if (!data.online) {
      statusEl.textContent = "🔴 Offline";
      statusEl.className = "status offline";
      motdEl.textContent = "";
      return;
    }

    statusEl.textContent =
      `🟢 Online — ${data.players.online}/${data.players.max} players`;
    statusEl.className = "status online";

    if (data.motd?.clean) {
      motdEl.textContent = data.motd.clean.join(" ");
    }

  } catch {
    statusEl.textContent = "⚠️ Unable to fetch server status";
    statusEl.className = "status offline";
  }
}

updateStatus();
setInterval(updateStatus, 60000); // refresh every 60s
