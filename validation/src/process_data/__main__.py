"""Allow `python -m process_data` from src/."""

import sys
from . import main

end_hour = None
if "--end-hour" in sys.argv:
    idx = sys.argv.index("--end-hour")
    end_hour = float(sys.argv[idx + 1])

main(end_hour=end_hour)
