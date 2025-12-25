const SERVER = "mc.cjones.dev";

const statusEl = document.getElementById("mc-status");
const playersEl = document.getElementById("mc-players");

async function updateStatus() {
  if (!statusEl) return;

  try {
    const res = await fetch(`https://api.mcstatus.io/v2/status/java/${SERVER}`);
    const data = await res.json();

    if (!data.online) {
      statusEl.textContent = "🔴 Offline";
      statusEl.className = "status offline";

      if (playersEl) {
        playersEl.innerHTML = "";
        playersEl.classList.add("empty");
        playersEl.innerHTML = "<li>No players online</li>";
      }
      return;
    }

    // Online
    statusEl.textContent =
      `🟢 Online — ${data.players.online}/${data.players.max} players`;
    statusEl.className = "status online";

    // Player list
    if (playersEl) {
      playersEl.innerHTML = "";

      if (data.players.sample && data.players.sample.length > 0) {
        data.players.sample.forEach(player => {
          const li = document.createElement("li");
          li.textContent = player.name;
          playersEl.appendChild(li);
        });
      } else {
        playersEl.classList.add("empty");
        playersEl.innerHTML = "<li>No players visible</li>";
      }
    }

  } catch (err) {
    statusEl.textContent = "⚠️ Unable to fetch server status";
    statusEl.className = "status offline";
  }
}

updateStatus();
setInterval(updateStatus, 60000);
