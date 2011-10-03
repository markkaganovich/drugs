from sqlalchemy import create_engine, Table, Column, Integer, String, MetaData, ForeignKey

engine = create_engine('sqlite:///:lcls:', echo=True)

metadata = MetaData()
lcls = Table('lcls', metadata, 
    Column('id', Integer, primary_key=True),
    Column('projectstage', String),
    Column('ethnicity', String),
    Column('status', String)
)

metadata.create_all(engine)


def lclimport(id, projectstage='', ethnicity='', status=''):
    conn = engine.connect()
    ins = lcls.insert().values(id=id, projectstage=projectstage,
ethnicity=ethnicity, status=status)
    conn.execute(ins)            
