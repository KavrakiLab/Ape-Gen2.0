#!/bin/bash

echo "🧹 Starting disk cleanup for GitHub Codespace / Ubuntu..."

# Show disk usage before cleanup
echo -e "\n💾 Disk usage BEFORE cleanup:"
df -h /

# Clean apt cache and remove unused packages
echo -e "\n🧼 Cleaning apt cache..."
sudo apt-get clean
sudo rm -rf /var/lib/apt/lists/*
sudo apt-get autoremove -y

# Clean pip cache
echo -e "\n🧼 Cleaning pip cache..."
rm -rf ~/.cache/pip

# Clean conda or mamba cache if available
echo -e "\n🧼 Cleaning conda/mamba cache (if available)..."
command -v mamba && mamba clean -a -y || true
command -v conda && conda clean -a -y || true

# Clean Docker data (if installed)
echo -e "\n🐳 Cleaning Docker images, containers, and volumes (if available)..."
command -v docker && docker system prune -a -f --volumes || true

# Remove VS Code remote cache
echo -e "\n🧽 Removing VS Code remote cache..."
rm -rf ~/.vscode-remote

# Remove Jupyter checkpoints
echo -e "\n📓 Removing Jupyter notebook checkpoint files..."
find ~ -type d -name ".ipynb_checkpoints" -exec rm -rf {} +

# Remove Python __pycache__ folders
echo -e "\n🐍 Removing Python __pycache__ folders..."
find ~ -type d -name "__pycache__" -exec rm -rf {} +

# Show disk usage after cleanup
echo -e "\n✅ Cleanup DONE! 💡 Disk usage AFTER cleanup:"
df -h /
