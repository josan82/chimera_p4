import cPickle, base64
try:
	from SimpleSession.versions.v65 import beginRestore,\
	    registerAfterModelsCB, reportRestoreError, checkVersion
except ImportError:
	from chimera import UserError
	raise UserError('Cannot open session that was saved in a'
	    ' newer version of Chimera; update your version')
checkVersion([1, 12, 41623])
import chimera
from chimera import replyobj
replyobj.status('Restoring session...', \
    blankAfter=0)
replyobj.status('Beginning session restore...', \
    blankAfter=0, secondary=True)
beginRestore()

def restoreCoreModels():
	from SimpleSession.versions.v65 import init, restoreViewer, \
	     restoreMolecules, restoreColors, restoreSurfaces, \
	     restoreVRML, restorePseudoBondGroups, restoreModelAssociations
	molInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVRFyaWJib25JbnNpZGVDb2xvcnECSwxOfYdVCWJhbGxTY2FsZXEDSwxHP9AAAAAAAAB9h1UJcG9pbnRTaXplcQRLDEc/8AAAAAAAAH2HVQVjb2xvcnEFSwxLAH1xBihLAV1xB0sBYUsCXXEISwJhSwNdcQlLA2FLBF1xCksEYUsFXXELSwVhSwZdcQxLBmFLB11xDUsHYUsIXXEOSwhhSwldcQ9LCWFLCl1xEEsKYUsLXXERSwthdYdVCnJpYmJvblR5cGVxEksMSwB9h1UKc3RpY2tTY2FsZXETSwxHP/AAAAAAAAB9h1UMbW1DSUZIZWFkZXJzcRRdcRUoTk5OTk5OTk5OTk5OZVUMYXJvbWF0aWNNb2RlcRZLDEsBfYdVCnZkd0RlbnNpdHlxF0sMR0AUAAAAAAAAfYdVBmhpZGRlbnEYSwyJfYdVDWFyb21hdGljQ29sb3JxGUsMTn2HVQ9yaWJib25TbW9vdGhpbmdxGksMSwB9h1UJYXV0b2NoYWlucRtLDIh9h1UKcGRiVmVyc2lvbnEcSwxLAH2HVQhvcHRpb25hbHEdfXEeVQhvcGVuZWRBc3EfiIlLDChVKi9ob21lL2luc2lsaWNoZW0vTGlnYW5kcy9SXzRrN2FfbGlnYW5kLnNkZnEgVQtNREwgTU9ML1NERnEhTol0cSJ9cSMoKFUoL2hvbWUvaW5zaWxpY2hlbS9MaWdhbmRzLzJhbTlfbGlnYW5kLnNkZnEkaCFOiXRxJV1xJksBYShVKi9ob21lL2luc2lsaWNoZW0vTGlnYW5kcy9SXzF4cTNfbGlnYW5kLnNkZnEnaCFOiXRxKF1xKUsKYShVKi9ob21lL2luc2lsaWNoZW0vTGlnYW5kcy9SXzFlM2tfbGlnYW5kLnNkZnEqaCFOiXRxK11xLEsJYShVKC9ob21lL2luc2lsaWNoZW0vTGlnYW5kcy8yaHZjX2xpZ2FuZC5zZGZxLWghTol0cS5dcS9LBmEoVSgvaG9tZS9pbnNpbGljaGVtL0xpZ2FuZHMvMmF4OV9saWdhbmQuc2RmcTBoIU6JdHExXXEySwVhKFUoL2hvbWUvaW5zaWxpY2hlbS9MaWdhbmRzLzJhbWJfbGlnYW5kLnNkZnEzaCFOiXRxNF1xNUsDYShVKC9ob21lL2luc2lsaWNoZW0vTGlnYW5kcy8yYXg2X2xpZ2FuZC5zZGZxNmghTol0cTddcThLBGEoVSgvaG9tZS9pbnNpbGljaGVtL0xpZ2FuZHMvM3Jsal9saWdhbmQuc2RmcTloIU6JdHE6XXE7SwdhKFUoL2hvbWUvaW5zaWxpY2hlbS9MaWdhbmRzLzFlM2dfbGlnYW5kLnNkZnE8aCFOiXRxPV1xPksAYShVKC9ob21lL2luc2lsaWNoZW0vTGlnYW5kcy8zcmxsX2xpZ2FuZC5zZGZxP2ghTol0cUBdcUFLCGEoVSgvaG9tZS9pbnNpbGljaGVtL0xpZ2FuZHMvMmFtYV9saWdhbmQuc2RmcUJoIU6JdHFDXXFESwJhdYeHc1UPbG93ZXJDYXNlQ2hhaW5zcUVLDIl9h1UJbGluZVdpZHRocUZLDEc/8AAAAAAAAH2HVQ9yZXNpZHVlTGFiZWxQb3NxR0sMSwB9h1UEbmFtZXFISwxYDwAAADJIVkNfTEdEX0FfMjIyNn1xSShYDwAAADFFM0tfUjE4X0FfMTAwMF1xSksJYVgPAAAAMUUzR19SMThfQV8xMDAwXXFLSwBhWAwAAAAzUkxMX1JMTF9BXzFdcUxLCGFYDwAAADFYUTNfUjE4X0FfMTAwMV1xTUsKYVgMAAAAMkFYOV9CSE1fQV8xXXFOSwVhWAwAAAAzUkxKX1JMSl9BXzFdcU9LB2FYDwAAADJBTUFfREhUX0FfMTAwMV1xUEsCYVgPAAAAMkFNQl8xN0hfQV8xMDAxXXFRSwNhWA8AAAAyQU05X1RFU19BXzEwMDBdcVJLAWFYDwAAADRLN0FfREhUX0FfMTAwMV1xU0sLYVgMAAAAMkFYNl9IRlRfQV8xXXFUSwRhdYdVD2Fyb21hdGljRGlzcGxheXFVSwyJfYdVD3JpYmJvblN0aWZmbmVzc3FWSwxHP+mZmZmZmZp9h1UKcGRiSGVhZGVyc3FXXXFYKH1xWX1xWn1xW31xXH1xXX1xXn1xX31xYH1xYX1xYn1xY31xZGVVA2lkc3FlSwxLCUsAhn1xZihLAEsAhl1xZ0sAYUsHSwCGXXFoSwdhSwNLAIZdcWlLA2FLCEsAhl1xaksIYUsGSwCGXXFrSwZhSwtLAIZdcWxLC2FLAksAhl1xbUsCYUsFSwCGXXFuSwVhSwpLAIZdcW9LCmFLAUsAhl1xcEsBYUsESwCGXXFxSwRhdYdVDnN1cmZhY2VPcGFjaXR5cXJLDEe/8AAAAAAAAH2HVRBhcm9tYXRpY0xpbmVUeXBlcXNLDEsCfYdVFHJpYmJvbkhpZGVzTWFpbmNoYWlucXRLDIh9h1UHZGlzcGxheXF1SwyIfYd1Lg=='))
	resInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQZpbnNlcnRxAksMVQEgfYdVC2ZpbGxEaXNwbGF5cQNLDIl9h1UEbmFtZXEESwxYAwAAAFVOS32HVQVjaGFpbnEFSwxYAQAAACB9h1UOcmliYm9uRHJhd01vZGVxBksMSwJ9h1UCc3NxB0sMiYmGfYdVCG1vbGVjdWxlcQhLDEsAfXEJKEsBTl1xCksBSwGGcQthhksCTl1xDEsCSwGGcQ1hhksDTl1xDksDSwGGcQ9hhksETl1xEEsESwGGcRFhhksFTl1xEksFSwGGcRNhhksGTl1xFEsGSwGGcRVhhksHTl1xFksHSwGGcRdhhksITl1xGEsISwGGcRlhhksJTl1xGksJSwGGcRthhksKTl1xHEsKSwGGcR1hhksLTl1xHksLSwGGcR9hhnWHVQtyaWJib25Db2xvcnEgSwxOfYdVBWxhYmVscSFLDFgAAAAAfYdVCmxhYmVsQ29sb3JxIksMTn2HVQhmaWxsTW9kZXEjSwxLAX2HVQVpc0hldHEkSwyJfYdVC2xhYmVsT2Zmc2V0cSVLDE59h1UIcG9zaXRpb25xJl1xJyhLAUsBhnEoSwFLAYZxKUsBSwGGcSpLAUsBhnErSwFLAYZxLEsBSwGGcS1LAUsBhnEuSwFLAYZxL0sBSwGGcTBLAUsBhnExSwFLAYZxMksBSwGGcTNlVQ1yaWJib25EaXNwbGF5cTRLDIl9h1UIb3B0aW9uYWxxNX1VBHNzSWRxNksMSv////99h3Uu'))
	atomInfo = cPickle.loads(base64.b64decode('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'))
	bondInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQVjb2xvcnECTS0BTn2HVQVhdG9tc3EDXXEEKF1xBShLGEsZZV1xBihLGEshZV1xByhLGUsaZV1xCChLGksbZV1xCShLGksrZV1xCihLG0scZV1xCyhLHEsdZV1xDChLHEshZV1xDShLHUseZV1xDihLHksfZV1xDyhLH0sgZV1xEChLH0slZV1xEShLIEshZV1xEihLIEsiZV1xEyhLIksjZV1xFChLI0skZV1xFShLJEslZV1xFihLJEsoZV1xFyhLJEspZV1xGChLJUsmZV1xGShLJksnZV1xGihLJ0soZV1xGyhLKEsqZV1xHChLKEssZV1xHShLLUsuZV1xHihLLUs3ZV1xHyhLLksvZV1xIChLL0swZV1xIShLL0sxZV1xIihLMUsyZV1xIyhLMkszZV1xJChLMks3ZV1xJShLM0s0ZV1xJihLNEs1ZV1xJyhLNUs2ZV1xKChLNUs7ZV1xKShLNks3ZV1xKihLNks4ZV1xKyhLN0tBZV1xLChLOEs5ZV1xLShLOUs6ZV1xLihLOks7ZV1xLyhLOks+ZV1xMChLOktAZV1xMShLO0s8ZV1xMihLPEs9ZV1xMyhLPUs+ZV1xNChLPks/ZV1xNShLQktDZV1xNihLQktMZV1xNyhLQ0tEZV1xOChLREtFZV1xOShLREtGZV1xOihLRktHZV1xOyhLR0tIZV1xPChLR0tMZV1xPShLSEtJZV1xPihLSUtKZV1xPyhLSktLZV1xQChLSktQZV1xQShLS0tMZV1xQihLS0tNZV1xQyhLTEtWZV1xRChLTUtOZV1xRShLTktPZV1xRihLT0tQZV1xRyhLT0tTZV1xSChLT0tVZV1xSShLUEtRZV1xSihLUUtSZV1xSyhLUktTZV1xTChLU0tUZV1xTShLV0tdZV1xTihLV0teZV1xTyhLWEteZV1xUChLWEtfZV1xUShLWUtaZV1xUihLWUtgZV1xUyhLWktfZV1xVChLW0tcZV1xVShLW0tgZV1xVihLW0thZV1xVyhLXEtiZV1xWChLXUthZV1xWShLXktsZV1xWihLX0thZV1xWyhLYEtkZV1xXChLYktjZV1xXShLY0tkZV1xXihLY0tnZV1xXyhLY0toZV1xYChLZEtlZV1xYShLZUtmZV1xYihLZktnZV1xYyhLZ0trZV1xZChLZ0ttZV1xZShLaEtpZV1xZihLaktrZV1xZyhLbktvZV1xaChLb0twZV1xaShLb0txZV1xaihLb0tyZV1xayhLcktzZV1xbChLckt0ZV1xbShLc0t6ZV1xbihLdEt1ZV1xbyhLdEt4ZV1xcChLdUt2ZV1xcShLdUt3ZV1xcihLeEt5ZV1xcyhLeUt6ZV1xdChLekt7ZV1xdShLe0t8ZV1xdihLfEt9ZV1xdyhLfEt+ZV1xeChLfkt/ZV1xeShLfkuAZV1xeihLfkuBZV1xeyhLgkuDZV1xfChLg0uEZV1xfShLhEuFZV1xfihLhEuGZV1xfyhLhEuHZV1xgChLh0uIZV1xgShLh0uJZV1xgihLiUuKZV1xgyhLikuLZV1xhChLikuMZV1xhShLi0uSZV1xhihLjEuNZV1xhyhLjUuOZV1xiChLjkuPZV1xiShLjkuSZV1xiihLj0uQZV1xiyhLj0uRZV1xjChLkkuTZV1xjShLk0uUZV1xjihLk0uVZV1xjyhLk0uWZV1xkChLl0uYZV1xkShLl0ujZV1xkihLmEuwZV1xkyhLmEuZZV1xlChLmUuaZV1xlShLmkueZV1xlihLmkukZV1xlyhLm0ukZV1xmChLnEukZV1xmShLnUukZV1xmihLnkujZV1xmyhLnkufZV1xnChLn0ugZV1xnShLoEuhZV1xnihLoEumZV1xnyhLoUuiZV1xoChLokujZV1xoShLpUumZV1xoihLpUurZV1xoyhLpkunZV1xpChLp0uvZV1xpShLqEurZV1xpihLqUurZV1xpyhLqkurZV1xqChLrEuvZV1xqShLrUuvZV1xqihLrkuvZV1xqyhLsUuzZV1xrChLsUu4ZV1xrShLsku6ZV1xrihLs0u1ZV1xryhLtEu6ZV1xsChLtUu3ZV1xsShLtUu7ZV1xsihLtku6ZV1xsyhLt0u5ZV1xtChLt0u6ZV1xtShLuEu5ZV1xtihLuEu9ZV1xtyhLu0u8ZV1xuChLvUu+ZV1xuShLvku/ZV1xuihLvkvAZV1xuyhLwEvBZV1xvChLwEvCZV1xvShLwEvDZV1xvihLw0vEZV1xvyhLxEvFZV1xwChLxUvGZV1xwShLxUvKZV1xwihLxkvHZV1xwyhLx0vIZV1xxChLyEvJZV1xxShLyEvLZV1xxihLyUvKZV1xxyhLy0vMZV1xyChLzUvrZV1xyShLzkvVZV1xyihLz0vWZV1xyyhL0EvjZV1xzChL0UvrZV1xzShL0kvsZV1xzihL00vsZV1xzyhL1EvsZV1x0ChL1UvlZV1x0ShL1kvmZV1x0ihL10vYZV1x0yhL10vdZV1x1ChL2EveZV1x1ShL2UvaZV1x1ihL2UvkZV1x1yhL2kvlZV1x2ChL20vcZV1x2ShL20vmZV1x2ihL3EvnZV1x2yhL3UvpZV1x3ChL3kvqZV1x3ShL30vkZV1x3ihL30voZV1x3yhL4EviZV1x4ChL4EvrZV1x4ShL4UvjZV1x4ihL4UvkZV1x4yhL4kvnZV1x5ChL40vrZV1x5ShL5UvoZV1x5ihL5kvpZV1x5yhL50vqZV1x6ChL6EvsZV1x6ShL6UvqZV1x6ihL7UvuZV1x6yhL7Uv2ZV1x7ChL7kvvZV1x7ShL70vwZV1x7ihL700AAWVdce8oS/BL8WVdcfAoS/FL8mVdcfEoS/FL9mVdcfIoS/JL82VdcfMoS/NL9GVdcfQoS/RL9WVdcfUoS/RL+mVdcfYoS/VL9mVdcfcoS/VL92VdcfgoS/dL+GVdcfkoS/hL+WVdcfooS/lL+mVdcfsoS/lL/WVdcfwoS/lL/mVdcf0oS/pL+2Vdcf4oS/tL/GVdcf8oS/xL/WVdcgABAAAoS/1L/2VdcgEBAAAoS/1NAQFlXXICAQAAKE0CAU0DAWVdcgMBAAAoTQIBTQsBZV1yBAEAAChNAwFNBAFlXXIFAQAAKE0EAU0FAWVdcgYBAAAoTQQBTRUBZV1yBwEAAChNBQFNBgFlXXIIAQAAKE0GAU0HAWVdcgkBAAAoTQYBTQsBZV1yCgEAAChNBwFNCAFlXXILAQAAKE0IAU0JAWVdcgwBAAAoTQkBTQoBZV1yDQEAAChNCQFNDwFlXXIOAQAAKE0KAU0LAWVdcg8BAAAoTQoBTQwBZV1yEAEAAChNDAFNDQFlXXIRAQAAKE0NAU0OAWVdchIBAAAoTQ4BTQ8BZV1yEwEAAChNDgFNEgFlXXIUAQAAKE0OAU0TAWVdchUBAAAoTQ8BTRABZV1yFgEAAChNEAFNEQFlXXIXAQAAKE0RAU0SAWVdchgBAAAoTRIBTRQBZV1yGQEAAChNEgFNFgFlXXIaAQAAKE0XAU0YAWVdchsBAAAoTRcBTSEBZV1yHAEAAChNGAFNGQFlXXIdAQAAKE0ZAU0aAWVdch4BAAAoTRkBTRsBZV1yHwEAAChNGwFNHAFlXXIgAQAAKE0cAU0dAWVdciEBAAAoTRwBTSEBZV1yIgEAAChNHQFNHgFlXXIjAQAAKE0eAU0fAWVdciQBAAAoTR8BTSABZV1yJQEAAChNHwFNJQFlXXImAQAAKE0gAU0hAWVdcicBAAAoTSABTSIBZV1yKAEAAChNIQFNKwFlXXIpAQAAKE0iAU0jAWVdcioBAAAoTSMBTSQBZV1yKwEAAChNJAFNJQFlXXIsAQAAKE0kAU0oAWVdci0BAAAoTSQBTSoBZV1yLgEAAChNJQFNJgFlXXIvAQAAKE0mAU0nAWVdcjABAAAoTScBTSgBZV1yMQEAAChNKAFNKQFlZVUFbGFiZWxyMgEAAE0tAVgAAAAAfYdVCGhhbGZib25kcjMBAABNLQGIfYdVBnJhZGl1c3I0AQAATS0BRz/JmZmgAAAAfYdVC2xhYmVsT2Zmc2V0cjUBAABNLQFOfYdVCGRyYXdNb2RlcjYBAABNLQFLAX2HVQhvcHRpb25hbHI3AQAAfVUHZGlzcGxheXI4AQAATS0BSwJ9h3Uu'))
	crdInfo = cPickle.loads(base64.b64decode('gAJ9cQEoSwB9cQIoSwBdcQMoRz+21JnDIh8GR0A8X/C2pbgtR0AUbwJMpn4Uh3EERz/qIh6AGUl/R0A7ZeDG2mvWR0AXwrSZ7r7lh3EFR0AB1hro0X1+R0A7Rj5gonasR0AV9Xvucbiqh3EGR0AIQETpM7fzR0A8feEmMkjGR0AVKhtzrSZ5h3EHR0ACoD4BmSjCR0A9j2aQrCbCR0AT3MecQ2aQh3EIR0AIhj0nshD3R0A+23i5kY9dR0AS+ow+9gZph3EJR0AAhDfiKurbR0A//9DLiVLlR0AUEO2DLYv2h3EKRz/u3UtDs9SVR0A/1aP/B7gsR0AQFlSVQ22Nh3ELRz/JOSKTp5iCR0A+niBbQE+qR0ARZNwGXeTYh3EMRz/sQWUr0kIaR0A9jxk/LiOmR0ATMJEdAfiDh3ENR7/0ZjFSdA8IR0A+k5SyI3GTR0AQt51ANN3oh3EOR7/+u2keosZ4R0A/u0g0osB3R0ARQ9iw1TkBh3EPR7/zQYrBddJkR0BAfxrOP+TdR0ASfOvsNiGwh3EQRz+3o6QzyR+kR0BAhm0zgBXVR0APH8t7HVoWh3ERR7/UQ1e4LRtnR0BAoPE4i6zcR0ADu45cA8bZh3ESR7/4qm82zUJyR0BBEGakouwbR0AEjSStloigh3ETR7//iXrpQbRER0BBEpc4MSAFR0AQA7ckt6Wkh3EUR7/xdwfDyROcR0BAoWMFyozsR0AYYVhbqY3Ih3EVR8ALrXc8AIfbR0BA/kmzKVRPR0AQsbUDWCtmh3EWR0AF0G1f9Vl/R0A6KNNAra4WR0AVIuzCHzcyh3EXR7/5qVGTGQvQR0BBrhmfLLL/R0ASXVoEfuyIh3EYZVUGYWN0aXZlcRlLAHVLAX1xGihLAF1xGyhHP97wdB18lRVHQEEWBqfgS0tHQA5jEZEXnO6HcRxHv8lB8mxthPxHQEGMFsAtJhFHQAZ4uxdu41yHcR1Hv/nWmuVxZTJHQEFeCumrrW1HQATzY/k+CIKHcR5HwAJP3y3/0PhHQEGo0T/eQWhHP/wj3bjGkY6HcR9HwAG6bcgKhM9HQEDZemPydzpHQAuBnD7Yi2SHcSBHv/cwWcKU9RZHQEBgJDsN0uJHQA+Owjg9JeWHcSFHwACyQIjtJuJHQD+zs5hgyWRHQBMgIPq8BsmHcSJHv/RP6SzlbppHQD5x33Thcq9HQBJNdOaFMcKHcSNHP7wd/gKZCyhHQD6qvtL1jgtHQBRvolcfM0GHcSRHP+kYHIJYK5RHQD/prdGpYRNHQBKRCiWxx0SHcSVHv4xeNiKFuQBHQEBfjnfQc95HQA1idT+Vn0qHcSZHQAE5nAVwYa9HQD+qwc6TjHpHQBBCscBjSUqHcSdHQAgL6u3jMVtHQD7H9Ayog55HQBPareQ1ocSHcShHQAKg2WaayfBHQD2MYcjNU/RHQBWk6goVFU+HcSlHP+7YbLeewOdHQD2Bq4/RM7VHQBLsWukDTTGHcSpHP984N7V5625HQDwqf1Y1hv9HQBSOlAATujSHcStHP/u7/8bQs6JHQDtNhf0MsK5HQBPHhwBN8Y6HcSxHQAdDXXozX6hHQDxDkW4YE0VHQBNmGHaCr8+HcS1HQA+eriLptX9HQDvay+G80rpHQBa5PAXnSzOHcS5HQAHl0xRWvzxHQD1idfg2er1HQBuX4/LOJlCHcS9HP9ShQdgqvUFHQEAkQBAI4fVHQAIV6C9sIHSHcTBlaBlLAHVLAn1xMShLAF1xMihHP8bl5BQaT/hHQEBZTl9me7RHQAS28VqNlSqHcTNHv99e+Gw01qJHQEDl2cMC3yNHP/zrRSHYPiSHcTRHv/sQ+dYhOVxHQEEq0zavFYlHQAOb5C6P1eSHcTVHwAP4KtbkXdlHQEF3fXWvxitHP/weMVJSZoyHcTZHwAAL1GaO2ptHQEEaXxGan6NHQA8jFqc0mk6HcTdHv/KkmDcouCpHQECTrITWmXNHQBJI8/UCutWHcThHv/2WpEl4l4BHQD/X1fkwTEdHQBHBo+U/0FSHcTlHv/HEgKu6WTpHQD6nSTm/myZHQBOU+eDtlxyHcTpHP9bmS5dMteBHQD68yn+qVrJHQBQQ261BXo2HcTtHP+8eEkin1HpHQEAKQxc8IhZHQBOPdo+EVk2HcTxHP82jXc4XpxhHQECHbdiJwjpHQBAiZr9oCXSHcT1HQAOGJxeARslHQD/uiBi/sDhHQBHZhJGXqMWHcT5HQAe6vzxEDAJHQD6+++m8uJpHQBUi0oM3DRmHcT9HQAP/xeMdachHQD2BSC2Ak9NHQBI6WgKOtiCHcUBHP/EVzjbaUQ9HQD2ykjKxI8ZHQBCHe9ChgoeHcUFHP9f/PuePuDZHQDxfwuM6KFNHQBGIUP//D/qHcUJHP/d5gf1VeyZHQDtn88QUvhdHQBLL8K9xtjCHcUNHQAQdxvNp3YNHQDxGh3PZbUVHQBWdUlyyGw6HcURHQA4ZS5OxGnNHQDunKzkalZxHQBVnnIylAhiHcUVHQArylEe7GFtHQD1QG8spETFHQAqi/47JqNWHcUZHP+6qQ2Afr2xHQEEyjI6WNxFHQBDOSUK7RqOHcUdlaBlLAHVLA31xSChLAF1xSShHP+lo/gLqq65HQDtO+ZvAGUxHQBcC1tZ/lc2HcUpHQAfOQGI6dNZHQDx+QKMiWOFHQBVEok8YZ0CHcUtHQAHnGKRy2xhHQD/S6PEZehRHQA+X6gQEi7SHcUxHQAi+gUmYi3pHQD7czte+kz5HQBMqdpcGlseHcU1HP8hNtfThm9ZHQD65qo50rqFHQBJ1HnM+5V6HcU5Hv/Pru2s17UpHQD6yCn6n66hHQBGT5NZPOzuHcU9HP7s+MTQe1HxHQDxTlsRhdTBHQBNPYMIIlSqHcVBHQAHfNE88T6xHQDtB5dKHLDlHQBYPhXACuaiHcVFHQAK1XCy1v85HQD2X0rfQgZFHQBP9/+g7JBqHcVJHP+5t9AImU3pHQEAA0t0cNxhHQBKl8+M/G82HcVNHP+vL34NXhwhHQD2c0r8XtQNHQBNHw3NYhPCHcVRHv/8TmrXUF8xHQD/MzdoVJppHQBHcGYGf5niHcVVHv/SvbpFfgNhHQECLzsjTaYBHQBMQJ3uWOSGHcVZHP7l7gA82MChHQECTgnFqBNxHQBDwzAZduceHcVdHv7PaZSTxaGpHQECqnu5CyX5HQAYqLlvnx4CHcVhHv/f3DiUc9rpHQEDrw2l4Q59HQARzneSK80uHcVlHv//nqinG3x5HQEEU/fYRVxJHQA+XxFHWw1yHcVpHv/YJ5Vaz825HQECtoIjwacJHQBkAgBlFeKqHcVtHwAaXBvwGKmRHQEC0I5Z5FwVHQBrEDr3/2hCHcVxHwAKedpS3IxxHQEJK5iQ6BL9HQAvgzNVJB8yHcV1Hv/l7HO6gbHJHQEHGF+eutQlHQBFEyt8213GHcV5HQAaJ7rEvFQZHQDofe31xZFRHQBXxYET9xN2HcV9HwAraLz96alJHQED26DhD2nRHQBBVJF+jaDeHcWBlaBlLAHVLBH1xYShLAF1xYihHP+94gnedMDpHQDsjM60DeRJHQA+c2JeOdc2HcWNHP+ARu4YQz31HQDwKb6jtaKdHQBN5hX6PCP2HcWRHv+vstgLfhwFHQDwjR9JOBN1HQBJaZO8o1N6HcWVHP+T2nqwvOoFHQDuCSiZrKzRHQBicujRXfoqHcWZHP/Ics8XAxtRHQD1epEddATRHQBNbE/qn2fKHcWdHP9SkevRqf6ZHQD55mO1/1otHQBHfplioGFuHcWhHQAODWYyR8JFHQD2UPu6jzB9HQBSBsX1QgJ+HcWlHQAp4ai8DzHtHQDyMjftzcOJHQBYvxj+It8GHcWpHQBEk65/RewRHQDz3KJRu1PdHQBn6ZBlgWHiHcWtHQArbUeUnqBFHQDtgzGEZBCpHQBTSgUKQfBqHcWxHQAeWprfFthdHQD7lmAHHHKBHQBQtvQB7yF2HcW1HQAFzzxJKqONHQD/2XKKJmgdHQBLDMgFzjpKHcW5HP+piYz7UC25HQD++EdqvHzxHQBGTRkoUUIqHcW9HP7KYkJplTsdHQEBvPiLOnYNHQBBEijllVCyHcXBHv/UvKtyHYo9HQEB7ztoT9C1HQBD7tggtVeyHcXFHv//X1bO2YuFHQEAENtwhXhhHQBLqBLCc0G6HcXJHwACD9EzZROhHQEEcBfSYRvBHQA7WaWFPdwaHcXNHv/8e0o7F8hNHQEE5hQ4mcntHQALgUI4F+VyHcXRHv/e2kHpedw1HQEG04+MENLdHQBJGIduRbmqHcXVHwAtTCtolX6dHQEEJH4oTT7xHQBCD4WoMnuiHcXZlaBlLAHVLBX1xdyhLAF1xeChHwAURlH3eg/5HQEIBfRZLSDxHP/xzAFNuCsiHcXlHv/zrKdvlncBHQEEwP6MxrbBHQAMieOtnlzSHcXpHwABB1720V5FHQEEdZ7sp3z9HQA8ToQJNh6aHcXtHwAsRI6XaIHlHQEEYPh2i2JZHQBB4Fs2KAQaHcXxHv/TZiIXfDAhHQEGm/Gk55bdHQBKFeNdpYQ+HcX1Hv/aKtZaz04JHQEBy90849vlHQBDv/0djkhSHcX5HwAEf3MeVb/dHQD/9VTYdju5HQBKHq8Hcf8+HcX9Hv4IhKUY5FsBHQEBmXdszH2dHQBBJyU3+Wg+HcYBHP+m2R4S+l9RHQD+q5TuIQcBHQBFI3n6W4CWHcYFHP9YTm12TmGlHQD5j1gT5JLNHQBIGAw/kDMyHcYJHQAGIi4noPQZHQD/Xd9SZDqtHQBGiXPzgFuqHcYNHQAiTjDIWxHdHQD7XGblnmu9HQBKiNHIGYwmHcYRHQAR9h7laFgZHQD2Kl5kor+lHQBNeS3OLjbiHcYVHQAv9wpSKyIxHQDyHaApRhaZHQBRD8G63IS2HcYZHQBKYf/luDeJHQDyri8l9xZNHQBIdxzCtnqSHcYdHQApdqe7ZNOxHQDtrt25pxepHQBckT5U3gySHcYhHP/M0Z16ucNpHQD1T5a5aS3BHQBMOVV4SiqOHcYlHP+IiqbCeDq5HQDv/J3g+ewBHQBPDdf/NGayHcYpHP/GtfaIWD0JHQDsdc0FCTgRHQBAWPD47uR6HcYtHP+MuOQ9Y4VJHQDuQYbeBDWRHQBjlLZ9sgLaHcYxHv+lVAujnfFZHQDwRgPAp8pBHQBI/47E0Y76HcY1laBlLAHVLBn1xjihLAF1xjyhHQA2TNLop+6xHQD55SWuOcmVHQBhkSg4K+d6HcZBHQA1KuRm6D+xHQD0v+wWnQ/tHQBoR+ydH0dWHcZFHQAUfLrd8q0hHQDxXWYdSUNpHQBiUbQxbY5aHcZJHP/prbOalxzVHQDzdkS43kNhHQBVI0Yg30jaHcZNHv+bNC+0nogJHQDx4MMyOvhRHQBWQESJmAYeHcZRHP9+Z7/OHOvhHQDu4j1sHy0VHQAyTB8hsev6HcZVHP+XiFvfVnjNHQDrDS5RpjrNHQBYjl0VHPqiHcZZHP/ryYF3zYaBHQD4xFjZleFFHQBOLOIPjpFSHcZdHP+bgSkIpJnpHQD66jkZ5E7ZHQBBVSiHYXzOHcZhHP+itZfHbSHdHQEAIETruR3ZHQA1rcIxJbRSHcZlHP/0gjTDgmQNHQEBpO3QNmiBHQBBl/H70Fo2HcZpHQAaT6JZlMIpHQEAokoGP0ABHQBOkIuJGyHmHcZtHQAYAdMLHjFRHQD79SewzzhJHQBU5ClR0HLWHcZxHP+DL5+mfy1dHQDwDgP27c8hHQBOzgfwKiq+HcZ1Hv/ja4vGzXbtHQEBdGvwHbt1HQAyqkOgxfOuHcZ5Hv9QRZttqLJhHQEBFCzQgHcVHQAbdWFR6d56HcZ9Hv9LQAi9IarxHQEBoOtTJOfJHP/ejgXULuEqHcaBHwAImj88yt7xHQEFtm/jq/KhHQAVhyH9iHxWHcaFHwAmoeWhFghpHQEECMt7bT4xHQBI2LQSXrOGHcaJHv/J22JGFijFHQEFnsWjnylpHQBKPwzLoGjeHcaNHwAAr6MGWVUZHQEEOjVo3s0FHQA4+WYu283WHcaRHP/bzJoA9ZMhHQD9MMkmzVsRHP+ckeMb3KjCHcaVHP+RK9IzIGxdHQEBo2C2xJ8xHv+Tjk0NvNuyHcaZHP/5jPkeAkKVHQEDFKhfsXO5HP++PF28IexCHcadHP+5rLIWxnjJHQEBLOOEoPV5HP+bVHWNbGnyHcahHQBI95QoTXWRHQDyvp0JyThpHQB0Wm2BBZpaHcallaBlLAHVLB31xqihLAF1xqyhHQAFpwoNLwC9HQD/phA1/NfpHQBKRNVu91F+HcaxHP+QBwwiuFUlHQDt1JtiGEO1HQBh6l6L3FGOHca1HQAd7PGUW25FHQD7Q09owQ7NHQBO7iPl4QlqHca5HP/BJa93KaaZHQDsnqpjq+l5HQA9540oOVdGHca9HQAL406JvF15HQD2RNCSmb45HQBQHZPQtQECHcbBHv+0GmVprVmNHQDwIBSnC6n9HQBJrFGSGDiCHcbFHP/CDxHS6+eBHQD1d8Nl2btlHQBMsU1oBrDaHcbJHP+lTgnmjce1HQD+/fTrakWtHQBGo2QenoRSHcbNHP9BjohpWlTpHQD5/d+UHbZhHQBH5YWqV84+HcbRHP9yUR4X3b+ZHQDwCktLqaDhHQBODwguSCTSHcbVHQAmfeTuc9CNHQDxtQj6n72RHQBVF7TXfXliHcbZHQA7yVFVo2bdHQDuDhgh4MQpHQBZMxHrqfuqHcbdHP6GzhNu/qOBHQEBz+7CSTqVHQBBY1p8EoAeHcbhHv/V6xO1m199HQEB/nJUF0ERHQBEAHlGp5L+HcblHwABfz8cQQ3tHQEAMRoKQmFxHQBLOK3e2K4yHcbpHwABdgWHo1qJHQEEkLKj0kLxHQA7TQy4nivmHcbtHwAs8gHSiJrhHQEESOLCIYhlHQBCdGzHP5OWHcbxHv/eKch9bwwhHQEGyV1UxlWhHQBKdK4yd/S+Hcb1Hv/1zfGTh6FRHQEE2L+ItNMFHQALR0IAv6syHcb5HwAN8HAAud5pHQEGymoAorDBHP/tKnV+SIYSHcb9HwAVDKb/uI8JHQEJTUQrzObFHP/brSRi4ZdyHccBHwABRe8qeymxHQELiWn2VobdHP//MdfZc4lqHccFHwALM/bPRMutHQEOJiw9lW1pHP/ocvgCv0CGHccJHwAqP3pram9lHQEOi+aeaxKNHP+W7wRlOipSHccNHwA+2k3HSTLpHQEMZr/S4TuBHP7iGTumnYkCHccRHwA0yEHeITkFHQEJ1/wUSNr9HP90elvlt16CHccVHwA0oXOt4t+tHQEROvZIP/ddHP9PBuM9UICiHccZHwA8oLA3KK81HQETaeD9AkhxHP5hhYrEOgICHccdlaBlLAHVLCH1xyChLAF1xyShHv/sUviTKykRHQEGRtQJ7lcNHQBLwm9xzyhmHccpHQBBHUVX56FhHQDuRIBLpsCNHQBU8IvzMJcCHcctHwA7hFwgoFY5HQETw6qjQ+SBHP8RhjkoiJPCHccxHwAJe7RW+cN9HQD/h6n2f7dVHQBD6XOPvQPCHcc1HwAyx+GdFTEhHQED5O7mE7qZHQA5a5etV+eqHcc5HP+3g5gO/Gl5HQDuSGTITpq5HQBnnQS4KJpaHcc9HP+ub3FHt8IJHQDsQOSodSGhHQBEFRg/WpZ+HcdBHv+tGTSVHaPBHQDwPQ50a6gpHQBUI/S4Brv+HcdFHQArPNtF38tRHQDxyfwPcvDVHQBRmTi7Kw+OHcdJHwAy/q5T62lxHQERlQQUR5FlHP9pUUJIju4CHcdNHwBMtgWZcoq1HQELUTpfe47tHwAAhlf17GQaHcdRHwBHpOzF15I5HQEIshKpEX+1Hv/t8U70/fCyHcdVHQABFmG4zb8hHQD/UYGOEQKlHQBEiTEPmc4mHcdZHQAcLGoEQbbpHQD7KExyV0JhHQBIgR+Wk6xeHcddHwAQkNPgKNUBHQEOQOPA3jb9HP/3uDNsZs6yHcdhHwAF5lsiH5qpHQELrLRggTZVHQAGXUezaAUCHcdlHwBGLjw6SeblHQENXwxhD7vNHv/MokounQSiHcdpHwA3yuiVHJNxHQEIJrIgNPTxHv+KqtRkjIUiHcdtHP8kmKg6UmUhHQD5grGEiQ9JHQBKad+l3gaGHcdxHv/ySK5y5nF1HQEFCmzGt24JHQAKEDO0dwoiHcd1Hv8SDqwK399BHQEBhKh9yknpHQBBLA2Dan26Hcd5HwAJM0rvz/hJHQEHEsFBtXwpHP/qgfQodpSGHcd9Hv/jkVdkEfGFHQEBtb0FvBhZHQBAhBPobzByHceBHP+SRl5JSVjJHQD+kFQsnTopHQBFbafBHEk6HceFHQAODmYmOO05HQD2GXtT4bAZHQBNfiFyGOH6HceJHwAn3XCmVvNRHQEO2SpDr55pHP+gDnVtIecSHceNHwAST5C3nn8RHQEJkGh3qiHFHP/aKjM64XieHceRHP/FNoufmL0NHQD1V+SJt0ExHQBOYipQK94OHceVHwA0kgipB0YhHQEMv7urGENpHv6vQ0MNAvsCHceZHwAqSsBsoTbJHQEKLe5EIKhZHP9Bclrq4/wCHcedHwAHYE2j6pYJHQEEXeBaQrxNHQA3BxHFq+dKHcehHP+DYAyU+z1xHQDwCfwf1RWZHQBTn51o+2B+HcellaBlLAHVLCX1x6ihLAF1x6yhHP8dEkcwzNz9HQDxG9XFkPz9HQBOmUigGKPCHcexHP/AajgYND7lHQDtBPr3smOtHQBBp9NjPu5eHce1HQAL3DoE3liZHQDs8+wDcVwFHQBKze30IhZ2Hce5HQAiplo9pO2lHQDx9ucoSlVJHQBOH08uZkVCHce9HQAL82scQ3PNHQD2a/1kUVoBHQBOd4GHlREOHcfBHQAgNvisxn4lHQD7o+9Ihn3JHQBRuHs5j6CKHcfFHP//v9vPaCEZHQEACf3ROItpHQBUy6vEb0M6HcfJHP+4xrNJ81NtHQD/np0tFzVdHQBDe6O4+pQOHcfNHP8ntvuztpMxHQD6jGqmB6rBHQBGYrLMqPyGHcfRHP+ys/3X8ZPZHQD2P5sDY+dJHQBLk4zpEPyKHcfVHv/PM6LvnFHZHQD6oGMM5H8tHQBDRkgyx+hSHcfZHv/+dGl758etHQD/CpwVfdD5HQBFcLHtRX5uHcfdHv/RrDG/Uqq5HQECBEOo6iYlHQBLQOWusFuSHcfhHP7eh64JK7jJHQECMmQsVUA5HQBBcOx1umVmHcflHv8b0y3gTXg1HQECvu8MZondHQAUAJFagAueHcfpHv/oSNzfry9NHQEDi0cmQMY1HQAQqFnIUnp6HcftHv//ZidTgLyFHQEEQzGDehZtHQA+Tx8wS1JWHcfxHv/VB1CAnZEVHQECePMu3IbBHQBivMEqmZwKHcf1Hv/Z/YwOqnG5HQEG8+tXSs/1HQBEsdwuxFhOHcf5HQAcvZsBJRp1HQDoifiZsHX5HQBPnDkG7pwiHcf9HwArLJKdxAvdHQEEPRxARgrhHQBCksp1rrhCHcgABAABlaBlLAHVLCn1yAQEAAChLAF1yAgEAAChHP8MKYoy+VkpHQDxLqckETO9HQBQhaQhr+cOHcgMBAABHP+w0fVuX7dBHQDtCWYdCnJBHQBCsrZIIT0KHcgQBAABHQAJPudr5RJ1HQDs6pxKuT8NHQBJfWppwYcaHcgUBAABHQAfuzby6nJhHQDx9MkSyGzxHQBNOWyQARbWHcgYBAABHQAKXEe3xpHNHQD2eseSfR6dHQBOdhe2OISOHcgcBAABHQAgmAXwWQA5HQD7iJE0LIRtHQBSM0Myv+lGHcggBAABHQADzhsBXM0RHQEAKesfkCIBHQBRePCBcRyqHcgkBAABHP+8o5BFwC65HQD/mHFSRDthHQBBug6ExHRGHcgoBAABHP8xBdgwutPRHQD6qQSnu5o1HQBFz237MYrWHcgsBAABHP+v8foJ12uNHQD2W+zuo4npHQBMDs+BFu6SHcgwBAABHv/PFQZb+x1NHQD6kK5jRBM9HQBCk+eBer0eHcg0BAABHv/5WHiL+5xZHQD/EhGGKMhRHQBFBOzHC9sKHcg4BAABHv/NYlLqPGuVHQECB+MjWzG9HQBKwdvGTExyHcg8BAABHP8HBsPP4ZwpHQECTHEZwDLVHQBAYXuoyNMeHchABAABHv7/uTTf0/ktHQEDJgcl15KBHQATYUaDzzLCHchEBAABHv/lPzYHWXMBHQED6y8Nc8DtHQARzRT5T0NeHchIBAABHv//WRvBpp/VHQEESnbG2m25HQA/o1xb/amKHchMBAABHv/GeAN/1UgtHQECjWWH43OFHQBiPBsMXOO2HchQBAABHv/mHTAGhVvhHQEG+Hea5qVVHQBIVeVL511yHchUBAABHQAbGIlVQShlHQDocNX/LJ1FHQBL4hEVjZfyHchYBAABHwArOKPnLOiRHQEDy7hJzJK9HQBBZEH2kWgqHchcBAABlaBlLAHVLC31yGAEAAChLAF1yGQEAAChHv/UiNNRf4aJHQECfu02TjsFHQBI56wIficaHchoBAABHwAKJfW0F4ahHQEDyS4jrjiFHQA1N+ItxQ+6HchsBAABHv/zu/skqznJHQEEk+fAZGyNHQALAejgJIB6HchwBAABHwAE+1iqycBBHQEGqLX6pIo1HP/x0D8cVP/yHch0BAABHv+pK7MQv5YlHQEC5twqVs71HP/lPEObAPmSHch4BAABHv8w9zrsSvSxHQEAycniXvehHQAPDIRTF30SHch8BAABHv/A/LDeRXTBHQD8lhaKHogFHQALGXOqiEMWHciABAABHv+rwTw83WadHQD4uW8LDykZHQAvsvD/TllqHciEBAABHP983eOfvDiJHQD5pp9oUwd1HQBB0j4Dq81eHciIBAABHP+WqRiiC1Q5HQD/Nqn3xurVHQBKaKLcV7kKHciMBAABHv6JtMstbsZBHQEBuRCG2r2ZHQA8jpl7S/UOHciQBAABHQAE33Mc2znFHQEACJSrjdXlHQBMwU9L6vBqHciUBAABHQAe0RT96ePtHQD7En4sSvTtHQBMeNKnWEi+HciYBAABHQAKuSDGg3VdHQD2T4q2YWBZHQBWfJQ4ZusiHcicBAABHP+sJiw9AJu9HQD19jwVXBTZHQBTmJT5ZMHSHcigBAABHP+BoCN03I05HQDwQXKsgeHVHQBQS+YBbJJuHcikBAABHP/tKKfEqjv5HQDuHh5yYmIpHQBERqG7pgtuHcioBAABHQAb/CMlnCx9HQDxXAjrFs2FHQBLRZ8VHwHGHcisBAABHQA1vKHx10tpHQDucItTEHzxHQBY8fDms1fKHciwBAABHQAU4aUMdI3hHQD2NSAjSJv5HQBuIdoUgs4SHci0BAABHP+siYEZe1GVHQEENNZNkdStHQA7v2f5FcWCHci4BAABlaBlLAHV1Lg=='))
	surfInfo = {'category': (0, None, {}), 'probeRadius': (0, None, {}), 'pointSize': (0, None, {}), 'name': [], 'density': (0, None, {}), 'colorMode': (0, None, {}), 'useLighting': (0, None, {}), 'transparencyBlendMode': (0, None, {}), 'molecule': [], 'smoothLines': (0, None, {}), 'lineWidth': (0, None, {}), 'allComponents': (0, None, {}), 'twoSidedLighting': (0, None, {}), 'customVisibility': [], 'drawMode': (0, None, {}), 'display': (0, None, {}), 'customColors': []}
	vrmlInfo = {'subid': (0, None, {}), 'display': (0, None, {}), 'id': (0, None, {}), 'vrmlString': [], 'name': (0, None, {})}
	colors = {u'': ((0.819608, 0, 0.309804), 1, u''), u'Ru': ((0.141176, 0.560784, 0.560784), 1, u'default'), u'gold': ((1, 0.843137, 0), 1, u'default'), u'Re': ((0.14902, 0.490196, 0.670588), 1, u'default'), u'Rf': ((0.8, 0, 0.34902), 1, u'default'), u'Ra': ((0, 0.490196, 0), 1, u'default'), u'Rb': ((0.439216, 0.180392, 0.690196), 1, u'default'), u'Rn': ((0.258824, 0.509804, 0.588235), 1, u'default'), u'Rh': ((0.0392157, 0.490196, 0.54902), 1, u'default'), u'Be': ((0.760784, 1, 0), 1, u'default'), u'Ba': ((0, 0.788235, 0), 1, u'default'), u'Bh': ((0.878431, 0, 0.219608), 1, u'default'), u'Bi': ((0.619608, 0.309804, 0.709804), 1, u'default'), u'Bk': ((0.541176, 0.309804, 0.890196), 1, u'default'), u'Br': ((0.65098, 0.160784, 0.160784), 1, u'default'), u'H': ((1, 1, 1), 1, u'default'), u'P': ((1, 0.501961, 0), 1, u'default'), u'Os': ((0.14902, 0.4, 0.588235), 1, u'default'), u'Ge': ((0.4, 0.560784, 0.560784), 1, u'default'), u'Gd': ((0.270588, 1, 0.780392), 1, u'default'), u'Ga': ((0.760784, 0.560784, 0.560784), 1, u'default'), u'Pr': ((0.85098, 1, 0.780392), 1, u'default'),
u'deep pink': ((1, 0.0784314, 0.576471), 1, u'default'), u'Pt': ((0.815686, 0.815686, 0.878431), 1, u'default'), u'Pu': ((0, 0.419608, 1), 1, u'default'), u'C': ((0.564706, 0.564706, 0.564706), 1, u'default'), u'Pb': ((0.341176, 0.34902, 0.380392), 1, u'default'), u'Pa': ((0, 0.631373, 1), 1, u'default'), u'Pd': ((0, 0.411765, 0.521569), 1, u'default'), u'Cd': ((1, 0.85098, 0.560784), 1, u'default'), u'Po': ((0.670588, 0.360784, 0), 1, u'default'), u'Pm': ((0.639216, 1, 0.780392), 1, u'default'), u'purple': ((0.627451, 0.12549, 0.941176), 1, u'default'), u'Hs': ((0.901961, 0, 0.180392), 1, u'default'), u'Ho': ((0, 1, 0.611765), 1, u'default'), u'Hf': ((0.301961, 0.760784, 1), 1, u'default'), u'Hg': ((0.721569, 0.721569, 0.815686), 1, u'default'), u'He': ((0.85098, 1, 1), 1, u'default'), u'Md': ((0.701961, 0.0509804, 0.65098), 1, u'default'), u'Mg': ((0.541176, 1, 0), 1, u'default'), u'K': ((0.560784, 0.25098, 0.831373), 1, u'default'), u'Mn': ((0.611765, 0.478431, 0.780392), 1, u'default'), u'O': ((1, 0.0509804, 0.0509804), 1, u'default'), u'Mt': ((0.921569, 0, 0.14902), 1, u'default'),
u'S': ((1, 1, 0.188235), 1, u'default'), u'W': ((0.129412, 0.580392, 0.839216), 1, u'default'), u'sky blue': ((0.529412, 0.807843, 0.921569), 1, u'default'), u'Zn': ((0.490196, 0.501961, 0.690196), 1, u'default'), u'plum': ((0.866667, 0.627451, 0.866667), 1, u'default'), u'Eu': ((0.380392, 1, 0.780392), 1, u'default'), u'Es': ((0.701961, 0.121569, 0.831373), 1, u'default'), u'Er': ((0, 0.901961, 0.458824), 1, u'default'), u'Ni': ((0.313725, 0.815686, 0.313725), 1, u'default'), u'No': ((0.741176, 0.0509804, 0.529412), 1, u'default'), u'Na': ((0.670588, 0.360784, 0.94902), 1, u'default'), u'Nb': ((0.45098, 0.760784, 0.788235), 1, u'default'), u'Nd': ((0.780392, 1, 0.780392), 1, u'default'), u'Ne': ((0.701961, 0.890196, 0.960784), 1, u'default'), u'dodger blue': ((0.117647, 0.564706, 1), 1, u'default'), u'Np': ((0, 0.501961, 1), 1, u'default'), u'Fr': ((0.258824, 0, 0.4), 1, u'default'), u'Fe': ((0.878431, 0.4, 0.2), 1, u'default'), u'Fm': ((0.701961, 0.121569, 0.729412), 1, u'default'), u'B': ((1, 0.709804, 0.709804), 1, u'default'), u'F': ((0.564706, 0.878431, 0.313725), 1, u'default'),
u'Sr': ((0, 1, 0), 1, u'default'), u'N': ((0.188235, 0.313725, 0.972549), 1, u'default'), u'Kr': ((0.360784, 0.721569, 0.819608), 1, u'default'), u'Si': ((0.941176, 0.784314, 0.627451), 1, u'default'), u'Sn': ((0.4, 0.501961, 0.501961), 1, u'default'), u'Sm': ((0.560784, 1, 0.780392), 1, u'default'), u'V': ((0.65098, 0.65098, 0.670588), 1, u'default'), u'Sc': ((0.901961, 0.901961, 0.901961), 1, u'default'), u'Sb': ((0.619608, 0.388235, 0.709804), 1, u'default'), u'Sg': ((0.85098, 0, 0.270588), 1, u'default'), u'Se': ((1, 0.631373, 0), 1, u'default'), u'Co': ((0.941176, 0.564706, 0.627451), 1, u'default'), u'Cm': ((0.470588, 0.360784, 0.890196), 1, u'default'), u'Cl': ((0.121569, 0.941176, 0.121569), 1, u'default'), u'Ca': ((0.239216, 1, 0), 1, u'default'), u'Cf': ((0.631373, 0.211765, 0.831373), 1, u'default'), u'Ce': ((1, 1, 0.780392), 1, u'default'), u'Xe': ((0.258824, 0.619608, 0.690196), 1, u'default'), u'Tm': ((0, 0.831373, 0.321569), 1, u'default'), u'light green': ((0.564706, 0.933333, 0.564706), 1, u'default'), u'Cs': ((0.341176, 0.0901961, 0.560784), 1, u'default'),
u'Cr': ((0.541176, 0.6, 0.780392), 1, u'default'), u'Cu': ((0.784314, 0.501961, 0.2), 1, u'default'), u'La': ((0.439216, 0.831373, 1), 1, u'default'), u'Li': ((0.8, 0.501961, 1), 1, u'default'), u'Tl': ((0.65098, 0.329412, 0.301961), 1, u'default'), u'Lu': ((0, 0.670588, 0.141176), 1, u'default'), u'Lr': ((0.780392, 0, 0.4), 1, u'default'), u'Th': ((0, 0.729412, 1), 1, u'default'), u'Ti': ((0.74902, 0.760784, 0.780392), 1, u'default'), u'tan': ((0.823529, 0.705882, 0.54902), 1, u'default'), u'Te': ((0.831373, 0.478431, 0), 1, u'default'), u'Tb': ((0.188235, 1, 0.780392), 1, u'default'), u'Tc': ((0.231373, 0.619608, 0.619608), 1, u'default'), u'Ta': ((0.301961, 0.65098, 1), 1, u'default'), u'Yb': ((0, 0.74902, 0.219608), 1, u'default'), u'Db': ((0.819608, 0, 0.309804), 1, u'default'), u'Dy': ((0.121569, 1, 0.780392), 1, u'default'), u'At': ((0.458824, 0.309804, 0.270588), 1, u'default'), u'I': ((0.580392, 0, 0.580392), 1, u'default'), u'salmon': ((0.980392, 0.501961, 0.447059), 1, u'default'), u'U': ((0, 0.560784, 1), 1, u'default'), u'Y': ((0.580392, 1, 1), 1, u'default'),
u'Ac': ((0.439216, 0.670588, 0.980392), 1, u'default'), u'Ag': ((0.752941, 0.752941, 0.752941), 1, u'default'), u'Ir': ((0.0901961, 0.329412, 0.529412), 1, u'default'), u'Am': ((0.329412, 0.360784, 0.94902), 1, u'default'), u'Al': ((0.74902, 0.65098, 0.65098), 1, u'default'), u'As': ((0.741176, 0.501961, 0.890196), 1, u'default'), u'Ar': ((0.501961, 0.819608, 0.890196), 1, u'default'), u'Au': ((1, 0.819608, 0.137255), 1, u'default'), u'Zr': ((0.580392, 0.878431, 0.878431), 1, u'default'), u'In': ((0.65098, 0.458824, 0.45098), 1, u'default'), u'Mo': ((0.329412, 0.709804, 0.709804), 1, u'default'), u'light gray': ((0.827451, 0.827451, 0.827451), 1, u'default')}
	materials = {u'': ((0.85, 0.85, 0.85), 30), u'default': ((0.85, 0.85, 0.85), 30)}
	pbInfo = {'category': [u'distance monitor'], 'bondInfo': [{'color': (0, None, {}), 'atoms': [], 'label': (0, None, {}), 'halfbond': (0, None, {}), 'labelColor': (0, None, {}), 'labelOffset': (0, None, {}), 'drawMode': (0, None, {}), 'display': (0, None, {})}], 'lineType': (1, 2, {}), 'color': (1, 16, {}), 'optional': {'fixedLabels': (True, False, (1, False, {}))}, 'display': (1, True, {}), 'showStubBonds': (1, False, {}), 'lineWidth': (1, 1, {}), 'stickScale': (1, 1, {}), 'id': [-2]}
	modelAssociations = {}
	colorInfo = (18, (u'deep pink', (1, 0.0784314, 0.576471, 1)), {(u'green', (0, 1, 0, 1)): [17], (u'Br', (0.65098, 0.160784, 0.160784, 1)): [15], (u'light green', (0.564706, 0.933333, 0.564706, 1)): [3], (u'dodger blue', (0.117647, 0.564706, 1, 1)): [8], (u'F', (0.564706, 0.878431, 0.313725, 1)): [13], (u'', (0.358033, 0.260402, 0.804281, 1)): [10], (u'N', (0.188235, 0.313725, 0.972549, 1)): [14], (u'', (0.123026, 0.337066, 0.083208, 1)): [11], (u'purple', (0.627451, 0.12549, 0.941176, 1)): [9], (u'gold', (1, 0.843137, 0, 1)): [7], (u'sky blue', (0.529412, 0.807843, 0.921569, 1)): [1], (u'tan', (0.823529, 0.705882, 0.54902, 1)): [0], (u'O', (1, 0.0509804, 0.0509804, 1)): [12], (u'plum', (0.866667, 0.627451, 0.866667, 1)): [2], (u'light gray', (0.827451, 0.827451, 0.827451, 1)): [5], (u'salmon', (0.980392, 0.501961, 0.447059, 1)): [4], (u'yellow', (1, 1, 0, 1)): [16]})
	viewerInfo = {'cameraAttrs': {'center': (0.10849997615814, 30.985500019073, 4.3005), 'fieldOfView': 25.350031531442, 'nearFar': (14.047373293742, -9.0454421618897), 'ortho': False, 'eyeSeparation': 50.8, 'focal': 4.3005}, 'viewerAttrs': {'silhouetteColor': None, 'clipping': False, 'showSilhouette': False, 'showShadows': False, 'viewSize': 10.357999968211, 'labelsOnTop': True, 'depthCueRange': (0.5, 1), 'silhouetteWidth': 2, 'singleLayerTransparency': True, 'shadowTextureSize': 2048, 'backgroundImage': [None, 1, 2, 1, 0, 0], 'backgroundGradient': [('Chimera default', [(1, 1, 1, 1), (0, 0, 1, 1)], 1), 1, 0, 0], 'depthCue': True, 'highlight': 0, 'scaleFactor': 1, 'angleDependentTransparency': True, 'backgroundMethod': 0}, 'viewerHL': 17, 'cameraMode': 'mono', 'detail': 1.5, 'viewerFog': None, 'viewerBG': None}

	replyobj.status("Initializing session restore...", blankAfter=0,
		secondary=True)
	from SimpleSession.versions.v65 import expandSummary
	init(dict(enumerate(expandSummary(colorInfo))))
	replyobj.status("Restoring colors...", blankAfter=0,
		secondary=True)
	restoreColors(colors, materials)
	replyobj.status("Restoring molecules...", blankAfter=0,
		secondary=True)
	restoreMolecules(molInfo, resInfo, atomInfo, bondInfo, crdInfo)
	replyobj.status("Restoring surfaces...", blankAfter=0,
		secondary=True)
	restoreSurfaces(surfInfo)
	replyobj.status("Restoring VRML models...", blankAfter=0,
		secondary=True)
	restoreVRML(vrmlInfo)
	replyobj.status("Restoring pseudobond groups...", blankAfter=0,
		secondary=True)
	restorePseudoBondGroups(pbInfo)
	replyobj.status("Restoring model associations...", blankAfter=0,
		secondary=True)
	restoreModelAssociations(modelAssociations)
	replyobj.status("Restoring camera...", blankAfter=0,
		secondary=True)
	restoreViewer(viewerInfo)

try:
	restoreCoreModels()
except:
	reportRestoreError("Error restoring core models")

	replyobj.status("Restoring extension info...", blankAfter=0,
		secondary=True)


try:
	import StructMeasure
	from StructMeasure.DistMonitor import restoreDistances
	registerAfterModelsCB(restoreDistances, 1)
except:
	reportRestoreError("Error restoring distances in session")


def restoreMidasBase():
	formattedPositions = {}
	import Midas
	Midas.restoreMidasBase(formattedPositions)
try:
	restoreMidasBase()
except:
	reportRestoreError('Error restoring Midas base state')


def restoreMidasText():
	from Midas import midas_text
	midas_text.aliases = {}
	midas_text.userSurfCategories = {}

try:
	restoreMidasText()
except:
	reportRestoreError('Error restoring Midas text state')


def restore_volume_data():
 volume_data_state = \
  {
   'class': 'Volume_Manager_State',
   'data_and_regions_state': [ ],
   'version': 2,
  }
 from VolumeViewer import session
 session.restore_volume_data_state(volume_data_state)

try:
  restore_volume_data()
except:
  reportRestoreError('Error restoring volume data')


def restore_cap_attributes():
 cap_attributes = \
  {
   'cap_attributes': [ ],
   'cap_color': None,
   'cap_offset': 0.01,
   'class': 'Caps_State',
   'default_cap_offset': 0.01,
   'mesh_style': False,
   'shown': True,
   'subdivision_factor': 1.0,
   'version': 1,
  }
 import SurfaceCap.session
 SurfaceCap.session.restore_cap_attributes(cap_attributes)
registerAfterModelsCB(restore_cap_attributes)

geomData = {'AxisManager': {}, 'CentroidManager': {}, 'PlaneManager': {}}

try:
	from StructMeasure.Geometry import geomManager
	geomManager._restoreSession(geomData)
except:
	reportRestoreError("Error restoring geometry objects in session")


def restoreSession_RibbonStyleEditor():
	import SimpleSession
	import RibbonStyleEditor
	userScalings = []
	userXSections = []
	userResidueClasses = []
	residueData = [(12, 'Chimera default', 'rounded', u'unknown'), (13, 'Chimera default', 'rounded', u'unknown'), (14, 'Chimera default', 'rounded', u'unknown'), (15, 'Chimera default', 'rounded', u'unknown'), (16, 'Chimera default', 'rounded', u'unknown'), (17, 'Chimera default', 'rounded', u'unknown'), (18, 'Chimera default', 'rounded', u'unknown'), (19, 'Chimera default', 'rounded', u'unknown'), (20, 'Chimera default', 'rounded', u'unknown'), (21, 'Chimera default', 'rounded', u'unknown'), (22, 'Chimera default', 'rounded', u'unknown'), (23, 'Chimera default', 'rounded', u'unknown')]
	flags = RibbonStyleEditor.NucleicDefault1
	SimpleSession.registerAfterModelsCB(RibbonStyleEditor.restoreState,
				(userScalings, userXSections,
				userResidueClasses, residueData, flags))
try:
	restoreSession_RibbonStyleEditor()
except:
	reportRestoreError("Error restoring RibbonStyleEditor state")

trPickle = 'gAJjQW5pbWF0ZS5UcmFuc2l0aW9ucwpUcmFuc2l0aW9ucwpxASmBcQJ9cQMoVQxjdXN0b21fc2NlbmVxBGNBbmltYXRlLlRyYW5zaXRpb24KVHJhbnNpdGlvbgpxBSmBcQZ9cQcoVQZmcmFtZXNxCEsBVQ1kaXNjcmV0ZUZyYW1lcQlLAVUKcHJvcGVydGllc3EKXXELVQNhbGxxDGFVBG5hbWVxDWgEVQRtb2RlcQ5VBmxpbmVhcnEPdWJVCGtleWZyYW1lcRBoBSmBcRF9cRIoaAhLFGgJSwFoCl1xE2gMYWgNaBBoDmgPdWJVBXNjZW5lcRRoBSmBcRV9cRYoaAhLAWgJSwFoCl1xF2gMYWgNaBRoDmgPdWJ1Yi4='
scPickle = 'gAJjQW5pbWF0ZS5TY2VuZXMKU2NlbmVzCnEBKYFxAn1xA1UHbWFwX2lkc3EEfXNiLg=='
kfPickle = 'gAJjQW5pbWF0ZS5LZXlmcmFtZXMKS2V5ZnJhbWVzCnEBKYFxAn1xA1UHZW50cmllc3EEXXEFc2Iu'
def restoreAnimation():
	'A method to unpickle and restore animation objects'
	# Scenes must be unpickled after restoring transitions, because each
	# scene links to a 'scene' transition. Likewise, keyframes must be 
	# unpickled after restoring scenes, because each keyframe links to a scene.
	# The unpickle process is left to the restore* functions, it's 
	# important that it doesn't happen prior to calling those functions.
	import SimpleSession
	from Animate.Session import restoreTransitions
	from Animate.Session import restoreScenes
	from Animate.Session import restoreKeyframes
	SimpleSession.registerAfterModelsCB(restoreTransitions, trPickle)
	SimpleSession.registerAfterModelsCB(restoreScenes, scPickle)
	SimpleSession.registerAfterModelsCB(restoreKeyframes, kfPickle)
try:
	restoreAnimation()
except:
	reportRestoreError('Error in Animate.Session')

def restoreLightController():
	import Lighting
	Lighting._setFromParams({'ratio': 1.25, 'brightness': 1.16, 'material': [30.0, (0.85, 0.85, 0.85), 1.0], 'back': [(0.3574067443365933, 0.6604015517481455, -0.6604015517481456), (1.0, 1.0, 1.0), 0.0], 'mode': 'two-point', 'key': [(-0.3574067443365933, 0.6604015517481455, 0.6604015517481456), (1.0, 1.0, 1.0), 1.0], 'contrast': 0.83, 'fill': [(0.2505628070857316, 0.2505628070857316, 0.9351131265310294), (1.0, 1.0, 1.0), 0.0]})
try:
	restoreLightController()
except:
	reportRestoreError("Error restoring lighting parameters")


def restoreRemainder():
	from SimpleSession.versions.v65 import restoreWindowSize, \
	     restoreOpenStates, restoreSelections, restoreFontInfo, \
	     restoreOpenModelsAttrs, restoreModelClip, restoreSilhouettes

	curSelIds =  []
	savedSels = []
	openModelsAttrs = { 'cofrMethod': 4 }
	windowSize = (512, 384)
	xformMap = {0: (((0, 0, 1), 0), (0, 0, 0), True), 1: (((0, 0, 1), 0), (0, 0, 0), True), 2: (((0, 0, 1), 0), (0, 0, 0), True), 3: (((0, 0, 1), 0), (0, 0, 0), True), 4: (((0, 0, 1), 0), (0, 0, 0), True), 5: (((0, 0, 1), 0), (0, 0, 0), True), 6: (((0, 0, 1), 0), (0, 0, 0), True), 7: (((0, 0, 1), 0), (0, 0, 0), True), 8: (((0, 0, 1), 0), (0, 0, 0), True), 9: (((0, 0, 1), 0), (0, 0, 0), True), 10: (((0, 0, 1), 0), (0, 0, 0), True), 11: (((0, 0, 1), 0), (0, 0, 0), True)}
	fontInfo = {'face': ('Sans Serif', 'Normal', 16)}
	clipPlaneInfo = {}
	silhouettes = {0: True, 1: True, 2: True, 3: True, 4: True, 5: True, 6: True, 7: True, 8: True, 9: True, 10: True, 11: True, 601: True}

	replyobj.status("Restoring window...", blankAfter=0,
		secondary=True)
	restoreWindowSize(windowSize)
	replyobj.status("Restoring open states...", blankAfter=0,
		secondary=True)
	restoreOpenStates(xformMap)
	replyobj.status("Restoring font info...", blankAfter=0,
		secondary=True)
	restoreFontInfo(fontInfo)
	replyobj.status("Restoring selections...", blankAfter=0,
		secondary=True)
	restoreSelections(curSelIds, savedSels)
	replyobj.status("Restoring openModel attributes...", blankAfter=0,
		secondary=True)
	restoreOpenModelsAttrs(openModelsAttrs)
	replyobj.status("Restoring model clipping...", blankAfter=0,
		secondary=True)
	restoreModelClip(clipPlaneInfo)
	replyobj.status("Restoring per-model silhouettes...", blankAfter=0,
		secondary=True)
	restoreSilhouettes(silhouettes)

	replyobj.status("Restoring remaining extension info...", blankAfter=0,
		secondary=True)
try:
	restoreRemainder()
except:
	reportRestoreError("Error restoring post-model state")
from SimpleSession.versions.v65 import makeAfterModelsCBs
makeAfterModelsCBs()

from SimpleSession.versions.v65 import endRestore
replyobj.status('Finishing restore...', blankAfter=0, secondary=True)
endRestore({})
replyobj.status('', secondary=True)
replyobj.status('Restore finished.')

