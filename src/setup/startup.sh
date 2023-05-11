#!/bin/bash

sudo apt update

# install chrome remote desktop 
curl -L -o chrome-remote-desktop_current_amd64.deb \
    https://dl.google.com/linux/direct/chrome-remote-desktop_current_amd64.deb
sudo DEBIAN_FRONTEND=noninteractive \
    apt-get install --assume-yes ./chrome-remote-desktop_current_amd64.deb

# xwindows system desktop environment - xfce in this case

sudo DEBIAN_FRONTEND=noninteractive \
    apt install --assume-yes xfce4 desktop-base dbus-x11 xscreensaver

# Configure Chrome Remote Desktop to use Xfce by default:
sudo bash -c 'echo "exec /etc/X11/Xsession /usr/bin/xfce4-session" > /etc/chrome-remote-desktop-session'

# Because there is no display connected to your instance, disable the display manager service on instance
sudo systemctl disable lightdm.service

# optional - installing chrome browser
curl -L -o google-chrome-stable_current_amd64.deb \
https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb
sudo apt install --assume-yes ./google-chrome-stable_current_amd64.deb


# create swap file
sudo dd if=/dev/zero of=/swapfile bs=1024 count=1M
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile

# install cmake and ccmake
# sudo apt-get -y install cmake 

# install gcc
sudo apt -y install build-essential

# install OpenSSL for cmake
sudo apt-get -y install libssl-dev

# install ncurses for ccmake
# sudo apt-get install libncurses5-dev

wget https://github.com/Kitware/CMake/releases/download/v3.26.3/cmake-3.26.3.tar.gz
tar -xvzf cmake-3.26.3.tar.gz
cd cmake-3.26.3
./configure
make
sudo make install

# sudo apt-get -y install cmake-curses-gui


#install/build itk

wget https://github.com/InsightSoftwareConsortium/ITK/releases/download/v5.3.0/InsightToolkit-5.3.0.tar.gz
tar -xvzf InsightToolkit-5.3.0.tar.gz
cd InsightToolkit-5.3.0
mkdir build
cd build
cmake ..
make
sudo make install




# install slider 3d
DEBIAN_FRONTEND=noninteractive sudo apt-get install libpulse-dev libnss3 libglu1-mesa 
sudo apt-get install --reinstall libxcb-xinerama0

cd folder
wget https://download.slicer.org/bitstream/63f5bee68939577d9867b4c7
tar -xvzf 63f5bee68939577d9867b4c7

# install python3

# install conda needed
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh
chmod +x Miniconda3-py310_23.3.1-0-Linux-x86_64.sh
echo "aef279d6baea7f67940f16aad17ebe5f6aac97487c7c03466ff01f4819e5a651 *Miniconda3-py310_23.3.1-0-Linux-x86_64.sh" | shasum --check
sudo ./Miniconda3-py310_23.3.1-0-Linux-x86_64.sh -b -p /opt/miniconda3

# needed to be run again for all users
/opt/miniconda3/bin/conda init

# restart shell
source ~/.bashrc 

# install snakemake
sudo /opt/miniconda3/bin/conda install -y -n base -c conda-forge mamba
sudo /opt/miniconda3/bin/mamba create -y -c conda-forge -c bioconda -n snakemake snakemake


# download repo

git clone https://github.com/mukundraj/broad-registration.git
cd broad-registration
mkdir build
cd build
cmake ..
cd ..

# install unzip
sudo apt install unzip


# set up user password
# sudo passwd $(whoami)
