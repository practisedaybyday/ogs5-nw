cmake_minimum_required( VERSION 2.8 )

set( TFILES
  TestDataPointInPolygon.gli
  )

UPDATE_MODEL_FILES ( 
  ${PROJECT_BINARY_DIR}/tests/data/GEO
  TFILES )
