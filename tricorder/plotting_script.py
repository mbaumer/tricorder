#!/u/ki/mbaumer/anaconda/bin/python
import subprocess

cmd_list = ["import results; res = results.Results('ZSPEC0.15_0.3nside1024nJack30','test4','redmagicHD',30); res.analyze_many_scales()",
            "import results; res = results.Results('ZSPEC0.3_0.45nside1024nJack30','test4','redmagicHD',30); res.analyze_many_scales()",
            "import results; res = results.Results('ZSPEC0.45_0.6nside1024nJack30','test4','redmagicHD',30); res.analyze_many_scales()",
            "import results; res = results.Results('DISTANCE0.15_0.3nside1024nJack30','test4','dark_matter',30); res.analyze_many_scales()",
            "import results; res = results.Results('DISTANCE0.3_0.45nside1024nJack30','test4','dark_matter',30); res.analyze_many_scales()",
            "import results; res = results.Results('DISTANCE0.45_0.6nside1024nJack30','test4','dark_matter',30); res.analyze_many_scales()",
            "import results; res = results.Results('DISTANCE0.6_0.75nside1024nJack30','test4','dark_matter',30); res.analyze_many_scales()",
            "import results; res = results.Results('DISTANCE0.75_0.9nside1024nJack30','test4','dark_matter',30); res.analyze_many_scales()"]

for cmd_str in cmd_list:
    subprocess.call(["bsub", "-W", "08:00", "-R", "rusage[mem=4000]",
                         "python", "-c", cmd_str])