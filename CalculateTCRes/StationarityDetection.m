function TC_KF_config = StationarityDetection( TC_KF_config,time )
%ÊµÏÖ¾²Ö¹¼ì²â
%   Detailed explanation goes here
if time>270099 && time<270345
    TC_KF_config.StationarityFlag=1;
else
    TC_KF_config.StationarityFlag=0;
end

end

