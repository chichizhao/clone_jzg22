% matlab 2022b
% author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
% script to plot the average non-synonymous ava
% plot non-synonymous avarage sample
data = [0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00904124453536231	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00911610417354073	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00891716095989163	0.00861876686266941	0.00861876686266941	0.00861876686266941	0.00893038348957199	0.00893038348957199	0.00893038348957199	0.00893038348957199	0.00893038348957199	0.00907808969902518	0.00881305823002550	0.00881305823002550	0.00903201832377742	0.00897742749868713	0.00918306672988546	0.00918306672988546	0.00918306672988546	0.00918306672988546	0.00918306672988546	0.00918306672988546	0.00918306672988546	0.00918306672988546	0.00918306672988546	0.00918306672988546	0.00918306672988546	0.00918306672988546	0.00918306672988546	0.00918306672988546	0.00934798118051420	0.00934798118051420	0.00934798118051420	0.00934798118051420	0.00934798118051420	0.00934798118051420	0.00934798118051420	0.00934798118051420	0.00934798118051420	0.00934798118051420	0.00934798118051420	0.00934798118051420	0.00934798118051420	0.00916107066751673	0.00916107066751673	0.00916107066751673	0.00916107066751673	0.00916107066751673	0.00916107066751673	0.00916107066751673	0.00960327058646084	0.00960327058646084	0.00960327058646084	0.00960327058646084	0.00960327058646084	0.00960327058646084	0.00960327058646084	0.00960327058646084	0.00968725479742917	0.00968725479742917	0.00968725479742917	0.00968725479742917	0.00968725479742917	0.00968725479742917	0.00968725479742917	0.00962684219517056	0.00962684219517056	0.00962684219517056	0.00962684219517056	0.00962684219517056	0.00962684219517056	0.00962684219517056	0.00962684219517056	0.00946109537993921	0.00946109537993921	0.00934803386226448	0.00934803386226448	0.00934803386226448	0.00934803386226448	0.00934803386226448	0.00934803386226448	0.00934803386226448	0.00934803386226448	0.00911342425265487	0.00911342425265487	0.00907032255292667	0.00907032255292667	0.00907032255292667	0.00907032255292667	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00893344497453507	0.00949196933069518	0.00949196933069518	0.00969348525801756	0.00969348525801756	0.00969348525801756	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0104236219876625	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.0102113420191609	0.00986277417243266	0.00986277417243266	0.00986277417243266	0.00918174415583314	0.00918174415583314	0.00918174415583314	0.00918174415583314	0.00918174415583314	0.00918174415583314	0.00918174415583314	0.00918174415583314	0.00939265377559190	0.00939265377559190	0.00939265377559190	0.00939265377559190	0.00939265377559190	0.00939265377559190	0.00939265377559190	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00971780109413369	0.00940623856062356	0.00940623856062356	0.00940623856062356	0.00940623856062356	0.00940623856062356	0.00940623856062356	0.00935663538602038	0.00935663538602038	0.00935663538602038	0.00883176638638228	0.00883176638638228	0.00883176638638228	0.00883176638638228	0.00883176638638228	0.00883176638638228	0.00883176638638228	0.00883176638638228	0.00883176638638228	0.00883176638638228	0.00883176638638228	0.00889072094004752	0.00820271846633076	0.00820271846633076	0.00820271846633076	0.00820271846633076	0.00694829651814876	0.00694829651814876	0.00694829651814876	0.00694829651814876	0.00694829651814876	0.00694829651814876	0.00694829651814876	0.00694829651814876	0.00631855569044503	0.00631855569044503	0.00631855569044503	0.00631855569044503	0.00631855569044503	0.00631855569044503	0.00631855569044503	0.00617119347411730	0.00617119347411730	0.00617119347411730	0.00617119347411730	0.00617119347411730	0.00617119347411730	0.00617119347411730	0.00626899042560725	0.00626899042560725	0.00626899042560725	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00610991977026051	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00593321634928228	0.00577870460638982	0.00577870460638982	0.00577870460638982	0.00577870460638982	0.00577870460638982	0.00577870460638982	0.00577870460638982	0.00577870460638982	0.00577870460638982	0.00577870460638982	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00689363146410106	0.00720124700688696	0.00720124700688696	0.00720124700688696	0.00720124700688696	0.00675751253759338	0.00675751253759338	0.00675751253759338	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00758980255934040	0.00770410072703345	0.00659083173770514	0.00659083173770514	0.00659083173770514	0.00659083173770514	0.00659083173770514	0.00659083173770514	0.00659083173770514	0.00659083173770514	0.00659083173770514	0.00659083173770514	0.00659083173770514	0.00609033123720464	0.00609033123720464	0.00609033123720464	0.00609033123720464	0.00609033123720464	0.00609033123720464	0.00609033123720464	0.00609033123720464	0.00609033123720464	0.00609033123720464	0.00609033123720464	0.00609033123720464	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00596506504654950	0.00399360270274673	0.00399360270274673	0.00532268297919543	0.00532268297919543	0.00532268297919543	0.00532268297919543	0.00532268297919543	0.00532268297919543	0.00532268297919543	0.00532268297919543	0.00532268297919543	0.00532268297919543	0.00532268297919543	0.00532268297919543	0.00532268297919543	0.00532268297919543	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00515534696848592	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00459531112619202	0.00444676686000069	0.00444676686000069	0.00444676686000069	0.00444676686000069	0.00444676686000069	0.00444676686000069	0.00444676686000069	0.00444676686000069	0.00444676686000069	0.00444676686000069	0.00444676686000069	0.00444676686000069	0.00444676686000069	0.00444676686000069	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00409662680397828	0.00388512088181246	0.00336101186713741	0.00336101186713741	0.00336101186713741	0.00336101186713741	0.00336101186713741	0.00336101186713741	0.00336101186713741	0.00194055732168287	0.00194055732168287	0.00194055732168287	0.00179638200449717	0.00179638200449717	0.00179638200449717	0.00179638200449717	0.00179638200449717	0.00179638200449717	0.00179638200449717	0.00179638200449717	0.00179638200449717	0.00179638200449717	0.0104452634185426	0.0104452634185426	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0190574008022317	0.0188758466991089	0.0188758466991089	0.0188758466991089	0.0188758466991089	0.0188758466991089	0.0188758466991089	0.0188758466991089	0.0165977878238278	0.0165977878238278	0.0165977878238278	0.0165977878238278	0.0165977878238278	0.0165977878238278	0.0165977878238278	0.0165977878238278	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0165182712843876	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.0263665108946972	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00757566122149456	0.00669382700632701	0.00669382700632701	0.00669382700632701	0.00669382700632701	0.00669382700632701	0.00669382700632701	0.00669382700632701	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.0306870548806033	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180	0.00486618004866180]
%set background color of figure
% region 1, from 1 to 1014
% region 2, from 1014 to 1114
% region 3, from 1114 to 1163
% region 4, from 1163 to 1194
% region 5, from 1194 to 1605
% here we use area function to plot the background color
% here we use 5 nice colors to plot the background color
color_set = [0.25 0.70 0.67; 1.00 0.50 0.31; 0.29 0.00 0.51; 0.42 0.56 0.14; 0.79 0.08 0.48];
color1 = color_set(1,:);
color2 = color_set(2,:);
color3 = color_set(3,:);
color4 = color_set(4,:);
color5 = color_set(5,:);
region1 = area([1 1014],[1 1],'FaceColor',color1,'EdgeColor',color1);
hold on
region2 = area([1014 1114],[1 1],'FaceColor',color2,'EdgeColor',color2);
hold on
region3 = area([1114 1163],[1 1],'FaceColor',color3,'EdgeColor',color3);
hold on
region4 = area([1163 1194],[1 1],'FaceColor',color4,'EdgeColor',color4);
hold on
region5 = area([1194 1605],[1 1],'FaceColor',color5,'EdgeColor',color5);
hold on
plot(data)
line([1014 1014],[0 1],'Color','k','LineStyle','--')
line([1114 1114],[0 1],'Color','k','LineStyle','--')
line([1163 1163],[0 1],'Color','k','LineStyle','--')
line([1194 1194],[0 1],'Color','k','LineStyle','--')
line([1605 1605],[0 1],'Color','k','LineStyle','--')
set(gca,'Xlabel',text('String','Time (years)','FontSize',18))
set(gca,'Ylabel',text('String','non-synonymous mutation number','FontSize',18))