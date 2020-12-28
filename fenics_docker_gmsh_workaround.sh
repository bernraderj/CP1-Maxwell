# Update the Ubuntu system.
sudo apt update
sudo apt dist-upgrade -y
# Install h5py, lxml and python3-pip.
sudo apt install -y python3-h5py python3-lxml python-h5py python-lxml python3-pip
# Update pip.
sudo pip install --upgrade pip
# Install pygmsh and meshio.
sudo -H pip  install -r requirements.txt
sudo -H pip3 install -r requirements.txt
# Install gmsh, latest version
sudo apt-get install -y wget libglu1 libxrender1  libxcursor-dev libxft-dev libxinerama-dev -y
sudo wget -nc  http://gmsh.info/bin/Linux/gmsh-4.6.0-Linux64-sdk.tgz -P /usr/local
sudo tar -xf /usr/local/gmsh-4.6.0-Linux64-sdk.tgz -C /usr/local/
sudo rm /usr/local/gmsh-4.6.0-Linux64-sdk.tgz
# Set variables.
echo "export PYTHONPATH=/usr/local/gmsh-4.6.0-Linux64-sdk/lib:$PYTHONPATH" >> ~/.profile
echo "export PATH=/usr/local/gmsh-4.6.0-Linux64-sdk/bin:$PATH" >> ~/.profile
echo "export PYTHONPATH=/usr/local/gmsh-4.6.0-Linux64-sdk/lib:$PYTHONPATH" >> ~/.bashrc
echo "export PATH=/usr/local/gmsh-4.6.0-Linux64-sdk/bin:$PATH" >> ~/.bashrc
# Replace 'webagg' with 'tkagg'.
sed -i 's/backend: webagg/backend: tkagg/g' ~/.config/matplotlib/matplotlibrc
# Then, make some alias commands.
echo "alias dir='ls -al'" >> ~/.bash_aliases
