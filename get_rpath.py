import site
print('-Wl,-rpath,.,-rpath,' + site.getusersitepackages() + '/mmoreseqs,-rpath,' + ',-rpath,'.join(site.getsitepackages()) + '/mmoreseqs')
