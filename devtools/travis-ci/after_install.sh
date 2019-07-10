# Temporarily change directory to $HOME to install software
pushd .
cd $HOME
git clone https://github.com/lpwgroup/respyte.git
export RESPYTEHOME=$HOME/respyte
cd $RESPYTEHOME
python setup.py install

# Restore original directory
popd
