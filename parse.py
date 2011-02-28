import pymongo.connection
import time

db = None

def createdb(file='CEUsnippet'):
    """This was your original function. I only removed the `eval`
    statements because they weren't necessary. When timed, this
    function took about 2 seconds when the database was empty."""
    global db
    f = open(file)
    lines = f.readlines()
    for line in lines:
        line = line.strip('\n')
        line = line.split('\t')
        if str(line[0]) == '#CHROM':
            cell_lines = line[9:]
            parameters = line[0:8]
        if '#' not in str(line[0]):
            record={}
            individuals={}
            for par in parameters:
                record[str(par)] = line[parameters.index(par)]
            for individual in cell_lines:
                individuals[str(individual)] = line[cell_lines.index(individual)+9]
                record.update(individuals)
            poslist = db['SNPs'+line[0]].distinct('POS')
            if line[1] in poslist:
                db['SNPs'+line[0]].update({'POS':line[1]}, {'$set':individuals})
            else:
                db['SNPs'+line[0]].insert(record)
    f.close()
    

def createdb2(filename='CEUSnippet'):
    """This is what I would do. You can simply do an upsert in this
    case: http://www.mongodb.org/display/DOCS/Updating#Updating-UpsertswithModifiers
    This function takes about 0.005 seconds, which is about 400x faster
    than the one above if my math serves me correctly."""
    global db
    f = open(filename)
    lines = f.readlines()
    keys = []
    for line in lines:
        line = line.strip('\n').split('\t')
        if line[0] == "#CHROM":
            keys = line
        if line[0].startswith('#'):
            continue
        record = {}
        for i, name in enumerate(keys):
            record[name.strip('#')] = line[i]
        db['SNPs'+record['CHROM']].update({'POS':record['POS']}, {'$set':record}, upsert=True)
    f.close()
        
        
    
if __name__ == "__main__":
    con = pymongo.connection.Connection()
    con.drop_database('test_db')
    db = pymongo.connection.Connection().test_db
    
    print "\nRunning your original function"
    start = time.time()
    try:
        createdb()
    except Exception:
        pass
    original_time = time.time() - start
    print "-" * 50
    print "Finished in %s seconds\n" % (original_time)
    
    print "Rebuilding the database\n"
    con.drop_database('test_db')
    db = pymongo.connection.Connection().test_db
    
    print "Running my version of the function"
    start = time.time()
    try:
        createdb2()
    except Exception:
        pass
    my_time = time.time() - start
    print "-" * 50
    print "Finished in %s seconds\n" % (my_time)
    
    print "=" * 50
    print "Updated function is about %sx faster\n" % int(original_time/my_time)