#!/u/ki/mbaumer/anaconda/bin/python
import subprocess

cmd_list = ["import results; res = results.Results('REDSHIFT0.15_0.3nside1024nJack30','test3','lss_sample',30); res.analyze_many_scales()",
            "import results; res = results.Results('REDSHIFT0.3_0.45nside1024nJack30','test3','lss_sample',30); res.analyze_many_scales()",
            "import results; res = results.Results('REDSHIFT0.45_0.6nside1024nJack30','test3','lss_sample',30); res.analyze_many_scales()",
            "import results; res = results.Results('REDSHIFT0.3_0.4nside1024nJack30','test3','lss_sample',30); res.analyze_many_scales()",
            "import results; res = results.Results('REDSHIFT0.6_0.7nside1024nJack30','test3','lss_sample',30); res.analyze_many_scales()",
            "import results; res = results.Results('REDSHIFT0.9_1.0nside1024nJack30','test3','lss_sample',30); res.analyze_many_scales()"
           ]

for cmd_str in cmd_list:
    subprocess.call(["bsub", "-W", "08:00", "-R", "rusage[mem=4000]",
                         "python", "-c", cmd_str])