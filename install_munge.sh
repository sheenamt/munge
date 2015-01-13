PREFIX=$1

python setup.py clean --all
rm -r $PREFIX/lib/python2.7/site-packages/munging
rm $PREFIX/bin/munge
python setup.py install --prefix=$PREFIX

