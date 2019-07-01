"""
Simple script retrieving annotations of given position.
"""

import os, sys, re, argparse, sqlite3, json
import Bio.Entrez, Bio.GenBank
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO

Bio.Entrez.email = 'your@email.address'

parser = argparse.ArgumentParser()
parser.add_argument('-g', default=None, metavar='string', help='Entrez ID')
parser.add_argument('-p', default=None, metavar='start[,stop]', help='Position')
parser.add_argument('-i', default=None, metavar='filename', help='Read GenBank file from file')
parser.add_argument('-o', default=None, metavar='filename', help='Output filenam (default:stdout)')
parser.add_argument('--db', default=None, metavar='filename', help='cache database')
parser.add_argument('--display', metavar='strings', default=None, nargs='+', help='display keywords such as gene, CDS (default show all features)')
parser.add_argument('--force', action='store_true', help='Force load from Entrez/local')
parser.add_argument('--verbose', action='store_true', help='verbosity')

indent = ' '  * 2
args = parser.parse_args()

forced = args.force
ostr = sys.stdout
if args.o is not None:
    ostr = open(args.o, 'w')

db = args.db
if db is None:
    dbdir = os.path.expanduser('~/.entrezcache')
    if os.path.exists(dbdir) is False:
        os.makedirs(dbdir)
    db = os.path.join(dbdir, 'genbank.cache.db')
entrez_id = args.g
pos = args.p
verbose = args.verbose
filename_input = args.i
display_key = args.display

if pos is not None and re.search('\\D', pos):
    start, stop = [int(x_.strip()) for x_ in re.split('\\D+', pos.strip())]
elif pos is not None:
    start = int(pos)
    stop = start
else:
    start = stop = 0

if verbose:
    if entrez_id is not None:
        sys.stderr.write('Genome ID:{}\n'.format(entrez_id))
    sys.stderr.write('Positoin:{}-{}\n'.format(start,stop))


with sqlite3.connect(db) as cnx:
    cur = cnx.cursor()
    cur.execute('SELECT name FROM sqlite_master WHERE name="genome"')
    if cur.fetchone() is None:
        cur.execute('CREATE TABLE genome (version NOT NULL PRIMARY KEY, accession, source, definition, size int8, contents blob)')
    cur.execute('SELECT name FROM sqlite_master WHERE name="feature"')
    if cur.fetchone() is None:
        cur.execute('CREATE TABLE feature (version NOT NULL, key, location, start int8, stop int8, qualifiers blob)')
        cur.execute('CREATE INDEX index_feature ON feature (version, start, stop)')
    cur.execute('SELECT version, size FROM genome WHERE version=? OR accession=?', (entrez_id, entrez_id))
    res = cur.fetchone()
    if res is None or forced:
        if filename_input is not None:
            contents = open(filename_input).read()
            rparser = Bio.GenBank.RecordParser()
            record = rparser.parse(StringIO(contents))
        else:
            handle = Bio.Entrez.efetch('nucleotide', id=entrez_id, rettype='gb', retmode='text')
            rparser = Bio.GenBank.RecordParser()
            record = rparser.parse(handle)
            contents = str(record)
            # with open('cache.txt', 'w') as fo:
            #     fo.write(contents.decode('utf8'))
            handle.close()
        base_length = 0 if not record.size.isdigit() else int(record.size)
        if forced:
            cur.execute('DELETE FROM genome WHERE version=?', (entrez_id,))
            cur.execute('DELETE FROM feature WHERE version=?', (entrez_id,))
        cur.execute('INSERT INTO genome VALUES(?, ?, ?, ?, ?, ?)', (record.version, ','.join(record.accession), record.source, record.definition, base_length, contents.encode('utf8')))
        cnx.commit()
        npat = re.compile('(\\d+)')
        for feature in record.features:
            key = feature.key
            location = feature.location
            m = npat.search(location)
            positions = []
            while m:
                positions.append(int(m.group(1)))
                m = npat.search(location, m.end())
            if len(positions) > 0:
                start = min(positions)
                stop = max(positions)
            else:
                start = stop = 0
                continue
            qualifiers = {}
            for qual in feature.qualifiers:
                if re.match('/(.+)=', qual.key):
                    label = qual.key[1:-1]
                    value = qual.value
                    qualifiers[key] = value.strip('"')
            cur.execute('INSERT INTO feature VALUES(?, ?, ?, ?, ?, ?)', (record.version, key, location, start, stop, json.dumps(qualifiers)))
        cnx.commit()
        version = record.version
    else:
        version = res[0]
    if start <= 0 and stop <= 0:
        start = 1
        stop = res[1]
    cur.execute('SELECT key, location, start, stop, qualifiers FROM feature WHERE version=? AND start <= ? AND stop >= ?', (version, stop, start))
    for key, location, start, stop, qualifiers in cur.fetchall():
        if display_key is None or key in display_key:
            output = '{} // {}\n'.format(key, location)
            qstr = ''
            for key, val in json.loads(qualifiers).items():
                qstr += '{}{} : {}\n'.format(indent, key, val)

            output += qstr
            ostr.write(output)

if args.o is not None:
    ostr.close()
