if [[ ! -z $1 ]]; then
PREFIX=$1
else
PREFIX=/home/genetics
fi

python setup.py clean --all
rm -r $PREFIX/lib/python2.7/site-packages/munging
rm $PREFIX/bin/munge
python setup.py install --prefix=$PREFIX

