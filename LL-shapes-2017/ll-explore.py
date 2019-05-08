from astropy.table import Table

df = Table(rows=intab[1:], names=intab[0]).to_pandas()
