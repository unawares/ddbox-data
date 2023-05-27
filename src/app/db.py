from app.configs import settings
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import NullPool

engine = create_engine(settings.DATABASE_URL, echo=False, future=True, pool_pre_ping=True, poolclass=NullPool)
Session = sessionmaker(bind=engine)
