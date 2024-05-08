import subprocess

# Install ffmpeg
subprocess.run(["sudo", "apt-get", "update"])
subprocess.run(["sudo", "apt-get", "install", "-y", "ffmpeg"])